from pathlib import Path
from os import system, path
import argparse

def get_stem(abs_path:Path) -> str:
    return path.splitext(path.basename(path.normpath(abs_path)))[0]

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Build script')
    parser.add_argument('src', type=Path, help='source file path')
    parser.add_argument('dest', type=Path, help='destination file path')
    parser.add_argument('--temp', type=bool, default=False,
                        help='automatically deletes the file if true')
    parser.add_argument('--compiler', default='g++',
                        help='compiler to use (default: g++)')
    parser.add_argument('--flags', default='-Wall -std=c++17',
                        help='compiler flags to use (default: -Wall -std=c++17)')

    args = parser.parse_args()
    abs_src = path.abspath(args.src)
    abs_dest = path.abspath(args.dest)


    system(f"{args.compiler} {args.flags} {abs_src} -o {abs_dest}/{args.src.stem}.o")
    # print(f"compiled in file {args.dest}")
    system(f'./{args.dest}/{args.src.stem}.o')

    if args.temp:
        system(f'rm {args.dest}/{get_stem(args.src)}.o')
