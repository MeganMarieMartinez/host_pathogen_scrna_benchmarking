set -eu

WORK=$HOME/dvsc_work

mkdir -p "$WORK/preproc" "$WORK/4" 
ln -sfn "$HOME/DVsc/src" "$WORK/src"  
ls -l "$WORK"