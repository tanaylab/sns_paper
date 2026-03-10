#!/usr/bin/env python
import os
import argparse
# Set backend
os.environ["KERAS_BACKEND"] = "torch"

from pathlib import Path
import numpy as np
import anndata
import crested
import keras
keras.mixed_precision.set_global_policy("mixed_bfloat16")
import zipfile
import tempfile
import glob
from crested.tl.zoo.utils._attention import MultiheadAttention  # force registration #added

def parse_args():
    parser = argparse.ArgumentParser(description='Train or finetune CRESTED model on genomic data')
    # Required arguments
    parser.add_argument('--bigwigs-folder', required=True, help='Path to folder containing bigwig files')
    parser.add_argument('--regions-file', required=True, help='Path to regions BED file')
    parser.add_argument('--output-dir', required=True, help='Output directory for results')
    parser.add_argument('--project-name', required=True, help='Project name for logging')
    
    # Training and model configuration
    parser.add_argument('--architecture', choices=['basenji', 'borzoi', 'chrombpnet', 'chrombpnet_decoupled', 
                                                'deeptopic_cnn', 'deeptopic_lstm', 'enformer', 'simple_convnet'], 
                        default='chrombpnet', help='Model architecture to use')
    parser.add_argument('--batch-size', type=int, default=128, help='Training batch size')
    parser.add_argument('--epochs', type=int, default=60, help='Number of training epochs')
    parser.add_argument('--logger', choices=['wandb', 'tensorboard'], default='wandb',
                       help='Logger to use for training')
    
    # Resume training
    parser.add_argument('--resume', action='store_true', help='Resume training from latest checkpoint')
    parser.add_argument('--checkpoint-path', help='Path to specific checkpoint to resume from (optional)')
    
    # Split strategy arguments
    parser.add_argument('--split-strategy', choices=['chr', 'region'], default='chr',
                       help='Strategy for train/val/test split: "chr" or "region"')
    parser.add_argument('--val-size', type=float, default=0.1,
                       help='Validation set size (only used if split_strategy="region")')
    parser.add_argument('--test-size', type=float, default=0.1,
                       help='Test set size (only used if split_strategy="region")')
    parser.add_argument('--val-chroms', nargs='+', default=['chr8', 'chr10'],
                       help='Validation chromosomes (only used if split_strategy="chr")')
    parser.add_argument('--test-chroms', nargs='+', default=['chr9', 'chr18'],
                       help='Test chromosomes (only used if split_strategy="chr")')
    
    # Borzoi finetuning specific arguments
    parser.add_argument('--finetune', action='store_true', help='Enable Borzoi finetuning mode')
    parser.add_argument('--borzoi-model', default='Borzoi_mouse_rep0', 
                        help='Borzoi model to use for finetuning (from crested.get_model)')
    parser.add_argument('--seq-len', type=int, default=2048, 
                        help='Sequence length for model input (must be divisible by 128 for Borzoi)')
    parser.add_argument('--initial-lr', type=float, default=1e-5, 
                        help='Initial learning rate for full dataset finetuning')
    parser.add_argument('--finetune-lr', type=float, default=5e-5, 
                        help='Learning rate for second-phase specific region finetuning')
    parser.add_argument('--finetune-epochs', type=int, default=5, 
                        help='Number of epochs for the second-phase finetuning')
    parser.add_argument('--gini-threshold', type=float, default=1.0, 
                        help='Threshold for filtering regions by specificity in second-phase finetuning')
    parser.add_argument('--two-phase', action='store_true', 
                        help='Enable two-phase training (first phase on all regions, second on specific regions)')
    parser.add_argument('--top-k-percent', type=float, default=0.03, 
                        help='Percent of top peaks to use for normalization')
    
    # Fixed paths that could be made configurable if needed
    parser.add_argument('--chromsizes-file', default='/home/aviezerl/proj/motif_reg/comparison/traj-for-crested/mm10.chrom.sizes',
                       help='Path to chromosome sizes file')
    parser.add_argument('--genome-file', default='/home/michalel/tanaylab/data/mm10/mm10.fa',
                       help='Path to genome FASTA file')
    
    return parser.parse_args()

def write_run_script(args, output_dir):
    """Write a shell script that can reproduce the run with the same parameters.
    
    Args:
        args: Parsed command line arguments
        output_dir: Path object for the output directory
    """
    script_path = output_dir / "reproduce_run.sh"
    
    # Get the path to the Python script itself
    python_script = os.path.abspath(__file__)
    
    # Start with the base command parts
    command_parts = [
        "#!/bin/bash",
        "",
        "# This script was automatically generated to reproduce the CRESTED run",
        "",
        f"python {python_script} \\"
    ]
    
    # Add required arguments (always present)
    required_args = [
        ("bigwigs-folder", args.bigwigs_folder),
        ("regions-file", args.regions_file),
        ("output-dir", args.output_dir),
        ("project-name", args.project_name),
        ("split-strategy", args.split_strategy),
        ("logger", args.logger),
        ("batch-size", args.batch_size),
        ("epochs", args.epochs),
        ("architecture", args.architecture)
    ]
    
    for arg_name, arg_value in required_args:
        command_parts.append(f"    --{arg_name} {arg_value} \\")
    
    # Add split strategy specific arguments
    if args.split_strategy == "region":
        command_parts.extend([
            f"    --val-size {args.val_size} \\",
            f"    --test-size {args.test_size} \\"
        ])
    else:  # chr strategy
        command_parts.extend([
            f"    --val-chroms {' '.join(args.val_chroms)} \\",
            f"    --test-chroms {' '.join(args.test_chroms)} \\"
        ])
    
    # Add finetuning specific arguments if enabled
    if args.finetune:
        command_parts.extend([
            f"    --finetune \\",
            f"    --borzoi-model {args.borzoi_model} \\",
            f"    --seq-len {args.seq_len} \\",
            f"    --initial-lr {args.initial_lr} \\",
            f"    --finetune-lr {args.finetune_lr} \\",
            f"    --top-k-percent {args.top_k_percent} \\"
        ])
        
        if args.two_phase:
            command_parts.extend([
                f"    --two-phase \\",
                f"    --finetune-epochs {args.finetune_epochs} \\",
                f"    --gini-threshold {args.gini_threshold} \\"
            ])
    
    # Add resume flag if it was used
    if args.resume:
        command_parts.append(f"    --resume \\")
        if args.checkpoint_path:
            command_parts.append(f"    --checkpoint-path {args.checkpoint_path} \\")
    
    # Add custom paths if they differ from defaults
    default_chromsizes = "/home/aviezerl/proj/motif_reg/comparison/traj-for-crested/mm10.chrom.sizes"
    default_genome = "/home/michalel/tanaylab/data/mm10/mm10.fa"
    
    if args.chromsizes_file != default_chromsizes:
        command_parts.append(f"    --chromsizes-file {args.chromsizes_file} \\")
    
    if args.genome_file != default_genome:
        command_parts.append(f"    --genome-file {args.genome_file} \\")
    
    # Remove the trailing backslash from the last line
    command_parts[-1] = command_parts[-1].rstrip(" \\")
    
    # Write the script
    with open(script_path, "w") as f:
        f.write("\n".join(command_parts))
    
    # Make the script executable
    os.chmod(script_path, 0o755)

def find_latest_checkpoint(output_dir, project_name, logger_type='wandb'):
    """Find the latest checkpoint file in the output directory.
    
    Args:
        output_dir: Path to the output directory
        project_name: Project name (for locating alternative checkpoint locations)
        logger_type: Type of logger used ('wandb' or 'tensorboard')
        
    Returns:
        Path to the latest checkpoint or None if no checkpoints found
    """
    # Try multiple potential checkpoint locations
    possible_checkpoint_dirs = [
        Path(output_dir) / logger_type / "checkpoints",  # Standard location
        Path(output_dir) / "checkpoints",                # Direct in output_dir
        Path(project_name) / Path(output_dir).name / "checkpoints"  # Project/run_name pattern
    ]
    
    # Try each location until we find checkpoints
    found_checkpoints = []
    for checkpoint_dir in possible_checkpoint_dirs:
        if not checkpoint_dir.exists():
            print(f"Checking {checkpoint_dir} - not found")
            continue
        
        # Look for checkpoint files (both numbered and best_model)
        checkpoint_files = list(checkpoint_dir.glob("*.keras"))
        
        if checkpoint_files:
            print(f"Found {len(checkpoint_files)} checkpoints in {checkpoint_dir}")
            found_checkpoints.extend(checkpoint_files)
    
    if not found_checkpoints:
        print("No checkpoint files found in any of the expected locations")
        return None
    
    # Filter out non-numbered checkpoints like best_model.keras
    numbered_checkpoints = [cp for cp in found_checkpoints if cp.stem.split(".")[0].isdigit()]
    
    if not numbered_checkpoints:
        # If we only have best_model.keras, return the first one found
        best_models = [cp for cp in found_checkpoints if cp.stem == "best_model"]
        if best_models:
            return best_models[0]
        return found_checkpoints[0]  # Return any checkpoint if no best_model
    
    # Sort numbered checkpoints by number
    sorted_checkpoints = sorted(numbered_checkpoints, 
                               key=lambda x: int(x.stem.split(".")[0]))
    
    # Return the latest checkpoint
    latest_checkpoint = sorted_checkpoints[-1]
    print(f"Found latest checkpoint: {latest_checkpoint}")
    return latest_checkpoint

def create_borzoi_scalar_model(seq_len, num_classes, pretrained_model_path=None, model_name=None):
    """Create a Borzoi model for scalar prediction finetuning.
    
    Args:
        seq_len: Sequence length (must be divisible by 128 for Borzoi)
        num_classes: Number of classes to predict
        pretrained_model_path: Optional path to pretrained model weights
        model_name: Name of the model (e.g., Borzoi_human_rep0, Borzoi_mouse_rep0)
        
    Returns:
        model_architecture: Keras model with Borzoi architecture modified for scalar prediction
    """
    # Determine number of pretrained model classes based on model type
    if model_name and "human" in model_name.lower():
        pretrained_num_classes = 7611  # Borzoi_human models have 7611 classes
    else:
        pretrained_num_classes = 2608  # Borzoi_mouse models have 2608 classes
    
    # Create default Borzoi architecture with specified input size
    # Target_length is set to input_len // 32 to match Borzoi's bin size
    base_model_architecture = crested.tl.zoo.borzoi(
        seq_len=seq_len, 
        target_length=seq_len//32, 
        num_classes=pretrained_num_classes
    )
    
    # Load pretrained weights if provided
    if pretrained_model_path:
        with zipfile.ZipFile(pretrained_model_path) as model_archive, tempfile.TemporaryDirectory() as tmpdir:
            model_weights_path = model_archive.extract('model.weights.h5', tmpdir)
            base_model_architecture.load_weights(model_weights_path)
    
    # Replace the head to predict scalar values per region
    # Get the output from the last layer of the base model before the head
    current = base_model_architecture.get_layer("final_conv_activation").output
    
    # Flatten the spatial dimension and add a new head
    current = keras.layers.Flatten()(current)
    current = keras.layers.Dense(
        num_classes, activation='softplus', name="dense_out"
    )(current)
    
    # Create a new model with the modified head
    model_architecture = keras.Model(
        inputs=base_model_architecture.inputs, 
        outputs=current, 
        name='Borzoi_scalar'
    )
    
    return model_architecture

def main():
    args = parse_args()
    
    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    run_name = output_dir.name
    
    # Check if resuming from checkpoint
    checkpoint_path = None
    if args.resume:
        if args.checkpoint_path:
            checkpoint_path = args.checkpoint_path
            if not os.path.exists(checkpoint_path):
                print(f"Warning: Specified checkpoint {checkpoint_path} not found.")
                checkpoint_path = None
        else:
            checkpoint_path = find_latest_checkpoint(output_dir, args.project_name, args.logger)
            
        if checkpoint_path:
            print(f"Resuming training from checkpoint: {checkpoint_path}")
        else:
            print("No valid checkpoint found. Starting training from scratch.")
            args.resume = False

    # Write the run script
    write_run_script(args, output_dir)
    
    # Skip preprocessing if resuming from checkpoint and data file exists
    adata_path = output_dir / "preprocessed_data.h5ad"
    if args.resume and adata_path.exists():
        print(f"Loading preprocessed data from {adata_path}")
        adata = anndata.read_h5ad(adata_path)
    else:
        # Import bigwigs
        adata = crested.import_bigwigs(
            bigwigs_folder=args.bigwigs_folder,
            regions_file=args.regions_file,
            target_region_width=200, #500
            chromsizes_file=args.chromsizes_file,
            target="mean" if not args.finetune else "mean", # Use "count" for ATAC-seq data
        )
        
        # Split data
        if args.split_strategy == "chr":
            crested.pp.train_val_test_split(
                adata,
                strategy="chr",
                val_chroms=args.val_chroms,
                test_chroms=args.test_chroms
            )
        else:
            crested.pp.train_val_test_split(
                adata,
                strategy="region",
                val_size=args.val_size,
                test_size=args.test_size,
                random_state=42
            )
        
        print(adata.var["split"].value_counts())
        
        # Preprocessing - adjust region width based on architecture and finetuning
        if args.finetune:
            target_width = args.seq_len
        else:
            # Default widths by architecture if not finetuning
            architecture_widths = {
                'basenji': 2114,
                'borzoi': 2048, # Must be divisible by 128
                'chrombpnet': 2114,
                'chrombpnet_decoupled': 2114,
                'deeptopic_cnn': 2114,
                'deeptopic_lstm': 2114,
                'enformer': 2114,
                'simple_convnet': 2114
            }
            target_width = architecture_widths.get(args.architecture, 2114)
        
        # Resize regions to target width
        crested.pp.change_regions_width(
            adata,
            target_width,
            chromsizes_file=args.chromsizes_file,
        )
        
        # Normalize peaks
        crested.pp.normalize_peaks(adata, top_k_percent=args.top_k_percent)
        
        # Save preprocessed data
        adata.write_h5ad(adata_path)
    
    # Setup data module
    datamodule = crested.tl.data.AnnDataModule(
        adata,
        genome=args.genome_file,
        batch_size=args.batch_size,
        max_stochastic_shift=3,
        always_reverse_complement=True,
    )
    
    # Model setup
    ##added
    if args.finetune and args.architecture.lower() == 'borzoi':
        if args.resume and checkpoint_path:
            # Load the model from checkpoint when resuming
            print(f"Loading model from checkpoint: {checkpoint_path}")
            #model_architecture = keras.models.load_model(checkpoint_path)
           
            from crested.tl.zoo.utils._attention import MultiheadAttention

            custom_objects = {
                "MultiheadAttention": MultiheadAttention,
            }

            model_architecture = keras.models.load_model(checkpoint_path, compile=True, custom_objects=custom_objects)
           
          #added  
        else:
            # Finetuning for Borzoi - create new model
            
            print(args.resume)
            
            print(checkpoint_path)
            
            print(f"Setting up Borzoi finetuning (seq_len={args.seq_len}, num_classes={len(list(adata.obs_names))})")
            
            # Get pretrained model
            model_file, _ = crested.get_model(args.borzoi_model)
            
            # Create Borzoi model with modified head for scalar prediction
            model_architecture = create_borzoi_scalar_model(
                seq_len=args.seq_len,
                num_classes=len(list(adata.obs_names)),
                pretrained_model_path=model_file,
                model_name=args.borzoi_model
            )
    else:
        if args.resume and checkpoint_path:
            # Load the model from checkpoint when resuming
            
            print(f"Loading model from checkpoint: {checkpoint_path}")
            model_architecture = keras.models.load_model(checkpoint_path)
        else:
            # Non-finetuning case - use regular model creation
            architecture_map = {
                'basenji': crested.tl.zoo.basenji,
                'borzoi': crested.tl.zoo.borzoi,
                'chrombpnet': crested.tl.zoo.chrombpnet,
                'chrombpnet_decoupled': crested.tl.zoo.chrombpnet_decoupled,
                'deeptopic_cnn': crested.tl.zoo.deeptopic_cnn,
                'deeptopic_lstm': crested.tl.zoo.deeptopic_lstm,
                'enformer': crested.tl.zoo.enformer,
                'simple_convnet': crested.tl.zoo.simple_convnet
            }
            
            # Get the selected architecture function
            architecture_func = architecture_map[args.architecture]
            
            # Initialize the model with the selected architecture
            model_architecture = architecture_func(
                seq_len=target_width,
                num_classes=len(list(adata.obs_names))
            )
    
    # Setup losses and metrics
    optimizer = keras.optimizers.Adam(learning_rate=args.initial_lr if args.finetune else 1e-3)
    loss = crested.tl.losses.CosineMSELogLoss(max_weight=100)
    metrics = [
        keras.metrics.MeanAbsoluteError(),
        keras.metrics.MeanSquaredError(),
        keras.metrics.CosineSimilarity(axis=1),
        crested.tl.metrics.PearsonCorrelation(),
        crested.tl.metrics.ConcordanceCorrelationCoefficient(),
        crested.tl.metrics.PearsonCorrelationLog(),
        crested.tl.metrics.ZeroPenaltyMetric(),
    ]
    
    config = crested.tl.TaskConfig(optimizer, loss, metrics)
    
    os.environ["KERAS_BACKEND"] = "torch"
    # First phase training
    print(f"Starting first phase training with learning rate: {optimizer.learning_rate.numpy()}")
    trainer = crested.tl.Crested(
        data=datamodule,
        model=model_architecture,
        config=config,
        project_name=args.project_name,
        run_name=run_name,
        logger=args.logger,
    )
    
    # Training
    # If resuming, we need to determine the starting epoch
    if args.resume and checkpoint_path:
        os.environ["KERAS_BACKEND"] = "tensorflow" #"tensorflow"
        # Extract the epoch number from the checkpoint filename if it's a numbered checkpoint
        checkpoint_name = os.path.basename(checkpoint_path)
        if checkpoint_name.split('.')[0].isdigit():
            start_epoch = int(checkpoint_name.split('.')[0])
            remaining_epochs = max(0, args.epochs - start_epoch)
            print(f"Resuming from epoch {start_epoch}. {remaining_epochs} epochs remaining.")
            
            if remaining_epochs > 0:
                trainer.fit(epochs=remaining_epochs, early_stopping=True, steps_per_epoch=4000)
            else:
                print("Training already completed. Skipping to evaluation.")
        else:
            # If using best_model or other non-numbered checkpoint, just run the full training
            # Early stopping will prevent unnecessary epochs
            print("Resuming from non-numbered checkpoint. Will use early stopping.")
            trainer.fit(epochs=args.epochs, early_stopping=True, steps_per_epoch=4000)
    else:        
        os.environ["KERAS_BACKEND"] = "tensorflow"#"tensorflow"
        # Regular training from scratch
        trainer.fit(epochs=args.epochs, early_stopping=True, steps_per_epoch=4000)
    
    # Testing
    trainer.test()    
   
    # Predictions
    trainer.predict(adata, model_name=run_name)
    
    # Save data after first phase
    adata.write_h5ad(f"{output_dir}/first_phase_predictions.h5ad")
    
    # Check if we should proceed with second phase
    # Either we're not resuming, or we're resuming and first phase is complete
    second_phase_file = f"{output_dir}/final_data_with_predictions.h5ad"
    should_run_second_phase = args.finetune and args.two_phase and (not args.resume or not os.path.exists(second_phase_file))
    
    # Second phase training (finetuning on specific regions) if enabled
    if should_run_second_phase:

        os.environ["KERAS_BACKEND"] = "tensorflow"#"torch" "tensorflow"
        print(f"Starting second phase finetuning on specific regions with threshold: {args.gini_threshold}")
        
        # Store original data for later predictions
        original_adata = adata.copy()
        
        # Filter regions based on specificity
        crested.pp.filter_regions_on_specificity(adata, gini_std_threshold=args.gini_threshold)
        print(f"After filtering, kept {adata.n_vars} regions")
        
        # Create new datamodule with filtered data
        datamodule = crested.tl.data.AnnDataModule(
            adata,
            genome=args.genome_file,
            batch_size=args.batch_size,
            max_stochastic_shift=3,
            always_reverse_complement=True,
        )        
        from keras.models import load_model#added
        from crested.tl.zoo.utils._attention import MultiheadAttention#added

        custom_objects = {
            "MultiheadAttention": MultiheadAttention,
        }#added

        

        # Load the best model from first phase
        best_model_path = find_latest_checkpoint(output_dir, args.project_name, args.logger)
        print(best_model_path) #added
        model = load_model(best_model_path, compile=True, custom_objects=custom_objects)#added
        #model = keras.models.load_model(best_model_path)
        
        # Update optimizer with new learning rate
        optimizer = keras.optimizers.Adam(learning_rate=args.finetune_lr)
        config = crested.tl.TaskConfig(optimizer, loss, metrics)
        
        os.environ["KERAS_BACKEND"] = "torch"# "tensorflow" torch
        # Setup new trainer for second phase
        ft_trainer = crested.tl.Crested(
            data=datamodule,
            model=model,
            config=config,
            project_name=args.project_name,
            run_name=f"{run_name}_phase2",
            logger=args.logger,
        )
        
        # Add this line before second phase training
        os.environ["KERAS_BACKEND"] = "tensorflow" # "tensorflow"
        # Train for fewer epochs in second phase
        ft_trainer.fit(epochs=args.finetune_epochs, early_stopping=True)
        
        # Test and predict with finetuned model
        ft_trainer.test()
        
        # Make predictions on all regions using the original data
        ft_trainer.predict(original_adata, model_name=f"{run_name}_phase2")
        
        # Update the main adata object with predictions from all regions
        adata = original_adata.copy()
    
    # Calculate and save correlations on test set
    # Skip if already done and we're resuming
    if not args.resume or not os.path.exists(f"{output_dir}/test_correlations.txt"):
        test_index = adata.var['split'] == 'test'
        correlations = []    
        
        # Determine which model name to use for correlations
        model_name = f"{run_name}_phase2" if (args.finetune and args.two_phase) else run_name
        
        with open(f"{output_dir}/test_correlations.txt", "w") as f:
            for traj in adata.obs_names:
                epi_x = adata.X[adata.obs_names == traj, test_index]
                epi_layer = adata.layers[model_name][adata.obs_names == traj, test_index]
                corr = np.corrcoef(epi_x.flatten(), epi_layer.flatten())[0][1]**2
                f.write(f"{traj} r^2: {corr}\n")
                correlations.append(corr)
            
            # Add mean correlation
            mean_corr = np.mean(correlations)
            f.write(f"\nMean r^2: {mean_corr}")
    
    # Save final data
    adata.write_h5ad(f"{output_dir}/final_data_with_predictions.h5ad")

if __name__ == "__main__":
    main()
