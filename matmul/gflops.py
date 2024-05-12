def gflops(n: int, time: int) -> float:
    return (2*n**3/(time*1e-9))/(1e9)  # 2n^3 flops / time (in ns) in GFLOPS


if __name__ == '__main__':
    n = 2048
    bench: dict = {
        'iterative': 634_936_604,
        'blocked': 19_087_074_146,
        'threaded_blocked': 11_113_214_750,
        'async_blocked': 10_820_800_292,
    }
    print(f'GFLOPS for {n}x{n} matrix multiplication')
    for name, time in bench.items():
        print(f'{name}: {gflops(n, time):.6f} GFLOPS')
