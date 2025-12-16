# コード19
load AF-P04637-F1-model_v4.pdb

# コード20
bg_color white

# コード21
select name , selection

# コード22
select AF_A179His, chain A and resi 179

# コード23
show spheres, AF_A179His

# コード24
run ../color_by_plddt.pml

# コード25
load pdb1gzh.ent

# コード26
select PDB_A179His, pdb1gzh and chain A and resi 179
show spheres, PDB_A179His

# コード27
center PDB_A179His

# コード28
select PDB_neighbors, pdb1gzh within 8 of PDB_A179His
show sticks, PDB_neighbors

# コード29
select ZnBinding, pdb1gzh and chain A and (resi 176 or resi 179 or resi 238 or resi 242)
hide sticks, PDB_neighbors
hide sphere, PDB_A179His
show sticks, ZnBinding

# コード30
color white, PDB_A179His and elem C

# コード31
label /pdb1gzh//A/HIS`179/CG, "His179"

# コード32
png pdb1ghk_His179.png, ray=1

# コード33
run ../hide_chainBCD.pml

# コード34
select PDB_chainA, pdb1gzh and chain A

# コード35
align AF-P04637-F1-model_v4, PDB_chainA

# コード36
select AF_ZnBinding, AF-P04637-F1-model_v4 and chain A and (resid 176 or resid 179 or resid 238 or resid 242
hide spheres, AF_A179His
show sticks, AF_ZnBinding

# コード37
run ../select_known_variants.pml

# コード38
color red, AF_A179His
show sphere, pathogenic

# コード39
hide sphere, pathogenic
show sphere, benign
show sphere AF_A179His

# コード50
load pdb7jjp.ent
bg_color white

# コード51
rotate X, 180

# コード52
select PDB_AtoF208Glu, resi 208 and (chain A or chain B or chain C or chain D or chain E or chain F)
show spheres, PDB_AtoF208Glu

# コード53
set cartoon_transparency, 0.5
set stick_transparency, 0.5

# コード54
hide everything, resn HOH

# コード55
rotate X, 90

# コード56
hide sticks, resn MC3

# コード57
select PDB_A208Glu, chain A and resid 208
center PDB_A208Glu
select neighbors, pdb7jjp within 8 of PDB_A208Glu
show sticks, neighbors

# コード58
hide sticks, neighbors
select Interactions, chain A and resi 205 or chain F and ( resi 76 or resi 73)
show sticks, PDB_A208Glu or Interactions

# コード59
set stick_transparency, 0

# コード60
label /pdb7jjp//A/GLU`208/CD, "Glu208"
label /pdb7jjp//F/ARG`76/CD, "Arg76"
label /pdb7jjp//A/ARG`205/CD, "Arg205"
label /pdb7jjp//F/SER`73/CA, "Ser73"

# コード61
set label_size, 20

# コード62
load AF-P48165-F1-model_v4.pdb

# コード63
select AF_Glu201, resid 201 and AF-P48165-F1-model_v4

# コード64
run ../color_by_plddt.pml

# コード65
hide cartoon, pdb7jjp and not chain A
hide labels
select PDBchainA, pdb7jjp and chain A
align AF-P48165-F1-model_v4, PDBchainA