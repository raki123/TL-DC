import sys

pw_tw_max_min_file = sys.argv[1]

result = [["path_width", "tree_width", "join_sum", "max_length", "min_length", "length_diff"]]

with open(pw_tw_max_min_file, "r") as stats:
    for line in stats:
        result.append([int(v) for v in line.split()])
        result[-1].append(result[-1][3] - result[-1][4])

for path in sys.argv[2:]:
    with open(path, "r") as solver_stats:
        result[0].append(path)
        for i, line in enumerate(solver_stats):
            line = line.split()
            if len(line) < 2:
                result[i + 1].append(-1)
            else:
                result[i + 1].append(int(line[1]))
with open("overall_stats.csv", "w") as out_file:
    for row in result:
        out_file.write(", ".join([str(val) for val in row]))
        out_file.write("\n")
