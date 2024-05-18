def gflops(n: int, time: int) -> float:
    return (2*n**3/(time*1e-9))/(1e9)  # 2n^3 flops / time (in ns) in GFLOPS


if __name__ == '__main__':
    n = 1024
    bench: dict = {
        'iterative': 84_804_104,
        'loop_reorder': 77_021_354,
        'blocked': 1_201_245_771,
        'threaded_gemm': 464_923_229,
        'neon': 478_612_646,
        'neon_threaded': 646_018_875,
    }
    print(f'GFLOPS for {n}x{n} matrix multiplication')
    for name, time in bench.items():
        print(f'{name}: {gflops(n, time):.6f} GFLOPS')
