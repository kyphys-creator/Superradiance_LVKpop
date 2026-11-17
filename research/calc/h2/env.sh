# Hoffman2(Slurm のクラスタ、hoffman2.oarc.ucla.edu)での Python の環境。ジョブの中で source する。
# 環境は ~/envs/srlvk(conda-forge の python 3.12.10・numpy 2.2.6・scipy 1.16.0・mpmath 1.3.0・sympy 1.14.0。Mac Studio と同じ版)。
# 公式の文書の作法(2026-11-17 に確認): モジュールは使うときに読み込む(~/.bashrc に書かない)、conda init はしない、
# PYTHONNOUSERSITE=1 で ~/.local の pip --user のパッケージが環境を上書きしないようにする。
# 1 つのジョブは 1 コアだけ使うので、BLAS のスレッドを 1 にする。
if ! type module >/dev/null 2>&1; then
  . /etc/profile.d/modules.sh 2>/dev/null || true
fi
module purge >/dev/null 2>&1
module load miniforge/26.5.3-0 >/dev/null 2>&1
conda activate "$HOME/envs/srlvk"
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export PYTHONWARNINGS=ignore
