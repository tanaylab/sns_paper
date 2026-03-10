mkdir -p ./mm10_zip

cd ./mm10 || exit 1

# mirror directory tree into ../mm10_zip
find . -type d -print0 | xargs -0 -I{} mkdir -p "../mm10_zip/{}"

# gzip each file into ../mm10_zip (keep originals in ./test_zip)
find . -type f ! -name '*.gz' -print0 | while IFS= read -r -d '' f; do
  gzip -c -9 -- "$f" > "../mm10_zip/$f.gz"
done