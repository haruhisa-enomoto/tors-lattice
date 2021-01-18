W = WeylGroup("A2", prefix='s')
s = W.simple_reflections()
W_pos = W.weak_lattice(facade = True)
j_irr = W_pos.join_irreducibles()
label = dict()
for A in W_pos.cover_relations_iterator():
    t = A[0] * A[1].inverse() # t: cover reflection ass. to v < w
    for x in j_irr:
        if W_pos.is_lequal(x,A[1]) and t in x.left_inversions_as_reflections():
            break
    label[(A[0],A[1])] = x
