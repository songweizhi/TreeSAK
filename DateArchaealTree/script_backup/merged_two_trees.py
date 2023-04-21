
t1 = ''
t2 = ''

leaves2n_a = {}
for n in t1.traverse():
    l = set(n.get_leaf_names())
    leaves2n_a[l] = n
leaves2n_b = {}
for n in t2.traverse():
    l = set(n.get_leaf_names())
    leaves2n_b[l] = n
    
n1_to_n2 = {}
for index, n1 in leaves2n_a.items():
    n2 = leaves2n_b[index]
    n1_to_n2[n1] = n2

## assume t1 has brancn length
## assume t2 has internal node name

merged_t = t1.copy()
for n, n1 in zip(merged_t.traverse(),
                 t1.traverse()):
    n.name = n1_to_n2[n1].name
merged_t.write(outfile='')