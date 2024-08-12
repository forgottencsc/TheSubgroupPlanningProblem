
# The Subgroup Planning Problem

## Dataset
```bash
curl -o - https://www.math.uwaterloo.ca/tsp/vlsi/vlsi_tsp.tgz | tar -zxvf -
```

## Dependencies

Install boost:

```bash
sudo apt install libboost-dev -y
```

## Compile

```bash
g++ tspp.cpp -o tspp -O3
```

## Run

```bash
./tspp <instance> <alg> <scale> <group by> <group size>
```

- `instance`: Path to instance
- `alg`: Which algorithm to run. Valid values: `1/2/3` for algorithms in this paper. `c` for comparison.
- `scale`: The multiple of edge weights inside each subgroup.
- `group by`: Method to generate subgroups. Valid values: `modulo|random|greedy(todo)`.
- `group size`: Maximum number of vertices for each group.

Output: A real number, denoting the weight of the tour obtained by the specified algorithm.

Example: 
```bash
./tspp vlsi/xql662.tsp c 3 modulo 10
```



