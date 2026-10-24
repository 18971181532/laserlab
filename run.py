"""LaserLab 统一引导脚本.

本机 Python 为嵌入式发行版（`python -m` 不含当前目录），故提供此
引导脚本，将项目根目录加入 `sys.path` 后再调用 CLI。

用法：

    python run.py cavity --length 0.2 --r1 0.99 --r2 0.98
    python run.py beam --length 0.2 --radius1 inf --radius2 0.5
    python run.py relax --pump-ratio 2.0
    python run.py modes --medium nd_yag_1064 --length 0.2
    python run.py serve --port 8779
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from laserlab.cli import main  # noqa: E402

if __name__ == "__main__":
    main()
