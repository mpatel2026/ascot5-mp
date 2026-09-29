#!/usr/bin/env python3
"""Summarize per-kernel register, stack and spill usage from a MINFO=1 build log.

Usage: tools/kernel_report.py variants/build.<tag>.log [filter]
       tools/kernel_report.py A.log B.log            (side-by-side comparison)
"""
import re
import sys

ENTRY = re.compile(r"Compiling entry function '(\w+)'")
FRAME = re.compile(r"(\d+) bytes stack frame, (\d+) bytes spill stores, "
                   r"(\d+) bytes spill loads")
REGS = re.compile(r"Used (\d+) registers")
# OpenACC kernel names look like _<len><file>_<function>_<line>_gpu
NAME = re.compile(r"_\d+(\w+?)_c_(\w+)_(\d+)_gpu$")


def parse(path):
    kernels = {}
    current = None
    frame = None
    with open(path) as f:
        for line in f:
            m = ENTRY.search(line)
            if m:
                current, frame = m.group(1), None
                continue
            m = FRAME.search(line)
            if m and current:
                frame = tuple(int(x) for x in m.groups())
                continue
            m = REGS.search(line)
            if m and current:
                name = current
                n = NAME.search(current)
                if n:
                    name = f"{n.group(2)}:{n.group(3)}"
                stack, st, ld = frame or (0, 0, 0)
                kernels[name] = (int(m.group(1)), stack, st, ld)
                current = None
    return {k: v for k, v in kernels.items() if "cub" not in k}


def main():
    logs = [a for a in sys.argv[1:] if a.endswith(".log")]
    filt = [a for a in sys.argv[1:] if not a.endswith(".log")]
    data = [parse(p) for p in logs]
    names = sorted(set().union(*data), key=lambda k: -data[0].get(k, (0,))[0])
    if filt:
        names = [n for n in names if any(f in n for f in filt)]
    head = "".join(f"{'regs stack spillST/LD':>28}" for _ in logs)
    print(f"{'kernel':50}{head}")
    for n in names:
        row = ""
        for d in data:
            if n in d:
                r, s, st, ld = d[n]
                row += f"{r:>8} {s:>6} {st:>6}/{ld:<6}"
            else:
                row += f"{'-':>28}"
        print(f"{n:50}{row}")


if __name__ == "__main__":
    main()
