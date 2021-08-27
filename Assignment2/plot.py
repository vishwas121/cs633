import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt

sns.set()

demo_input_format = pd.DataFrame.from_dict({
    "D": [],
    "P": [],
    "ppn": [],
    "mode": [],  # 1 --> optimized, 0 --> standard
    "time": [],
})

files = ['Bcast', 'Reduce', 'Gather', 'Alltoallv']
for file in files:
    f1 = open(f'data_{file}.txt', "r")
    data = f1.readlines()
    line = 0

    for execution in range(10):
        for P in [4, 16]:
            for ppn in [1, 8]:
                for D in [16, 256, 2048]:
                    nums = data[line].strip().split(" ")
                    if file is "Alltoallv" and D is 2048:
                        D /= 2
                    demo_input_format = demo_input_format.append({
                        "D": D, "P": P, "ppn": ppn, "mode": 0, "time": float(nums[0])
                    }, ignore_index=True)
                    demo_input_format = demo_input_format.append({
                        "D": D, "P": P, "ppn": ppn, "mode": 1, "time": float(nums[1])
                    }, ignore_index=True)
                    line += 1

    demo_input_format["(P, ppn)"] = list(map(lambda x, y: ("(" + x + ", " + y + ")"), map(str, demo_input_format["P"]), map(str, demo_input_format["ppn"])))


    g = sns.catplot(x="(P, ppn)", y="time", data=demo_input_format,
                kind="bar", col="D", hue="mode", sharey="col")
    g.tight_layout()
    plt.savefig(f'plot_{file}.jpg')