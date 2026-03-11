"""
Distributed computing utilities for compute_dependency_matrix.

This module contains distributed/multi-GPU coordination functions extracted from
compute_dependency_matrix.py for better code organization.
"""

import os
from typing import Tuple

import torch
import torch.distributed as dist


def setup_distributed() -> Tuple[int, int, int]:
    """
    Initialize distributed training if available.

    Returns:
        Tuple of (rank, world_size, local_rank)
    """
    if 'RANK' in os.environ and 'WORLD_SIZE' in os.environ:
        rank = int(os.environ['RANK'])
        world_size = int(os.environ['WORLD_SIZE'])
        local_rank = int(os.environ.get('LOCAL_RANK', 0))

        dist.init_process_group(backend='nccl')
        torch.cuda.set_device(local_rank)

        return rank, world_size, local_rank
    else:
        return 0, 1, 0


def cleanup_distributed() -> None:
    """Clean up distributed resources."""
    if dist.is_initialized():
        dist.destroy_process_group()


def print_rank0(msg: str, rank: int = 0) -> None:
    """Print only on rank 0."""
    if rank == 0:
        print(msg)


def synchronize_if_distributed() -> None:
    """Synchronize all processes if running in distributed mode."""
    if dist.is_initialized():
        dist.barrier()


def is_main_process(rank: int = 0) -> bool:
    """Check if this is the main process (rank 0)."""
    return rank == 0


def get_device_for_rank(local_rank: int) -> torch.device:
    """Get the appropriate CUDA device for this rank."""
    if torch.cuda.is_available():
        return torch.device(f'cuda:{local_rank}')
    return torch.device('cpu')


def all_gather_object(obj, world_size: int):
    """
    Gather objects from all ranks.

    Args:
        obj: Object to gather
        world_size: Number of processes

    Returns:
        List of objects from all ranks (on rank 0) or None (on other ranks)
    """
    if not dist.is_initialized() or world_size == 1:
        return [obj]

    gathered = [None] * world_size
    dist.all_gather_object(gathered, obj)
    return gathered


def broadcast_object(obj, src: int = 0):
    """
    Broadcast object from source rank to all ranks.

    Args:
        obj: Object to broadcast (only meaningful on src rank)
        src: Source rank

    Returns:
        Broadcasted object
    """
    if not dist.is_initialized():
        return obj

    object_list = [obj]
    dist.broadcast_object_list(object_list, src=src)
    return object_list[0]
