#export DYLD_LIBRARY_PATH=$CONDA_PREFIX/lib

python generate_sheared_structures.py
for i in {0..20};do
    mkdir shear-$i
    mv structure_$i shear-$i/unitcell
done

# Download machine learning potential
if [ ! -f mlp.lammps ]; then
    wget http://cms.mtl.kyoto-u.ac.jp/seko/mlp-repository/_downloads/02768f97a219d19cdbbf16d4f0736e2e/mlp.lammps
fi

# Relax unit cell
for i in {0..20};do
    cd shear-$i
    ~/code/lammps/src/lmp_conda -in ../in.relax
    cd ..
done

# Create supercells
for i in {0..20};do
    cd shear-$i
    phonopy -c relaxed_unitcell --lammps -d --dim 4 4 3 
    cd ..
done

# Run force calculations
for i in {0..20};do
    cd shear-$i
    for num in `ls supercell-*|sed s/supercell-//g`;do
        sed s/number/$num/ ../in.forces | ~/code/lammps/src/lmp_conda
    done
    cd ..
done

# Run phonono calculations
for i in {0..20};do
    cd shear-$i
    phonopy -f forces-*
    phonopy-load --band 0 0 1/2  1/3 1/3 1/2  1/3 1/3 0  0 0 0  1/2 0 0  1/2 0 1/2  0 0 1/2
    cd ..
done

phonopy-bandplot shear-{0,4,8,12,16,20}/band.yaml