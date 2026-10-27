"""Validate C++ against the original, frozen Python/JAX checkpoint evaluator."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time

os.environ.setdefault('JAX_PLATFORMS','cpu')
import jax
jax.config.update('jax_enable_x64',True)
import jax.numpy as jnp
import numpy as np


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workspace',type=Path,required=True,help='qsmhd-pinn repository root')
    parser.add_argument('--executable',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    root=args.workspace.resolve()/'free_boundary_equilibrium'
    source=root/'experiments/iter_taper15_expanded_20261026/source'
    snapshot=root/'analysis/expanded_current_47892285/step_10000'
    sys.path.insert(0,str(source))
    from train_expanded import ExpandedProblem
    from PINN_taper import taper, FB
    from PINN_free_boundary import FLUX, ICOIL
    saved=np.load(snapshot/'checkpoint.npz',allow_pickle=False)
    geometry=dict(np.load(source/'original_geometry.npz',allow_pickle=False))
    started=time.perf_counter()
    problem=ExpandedProblem(geometry,points=32)
    problem.mask=saved['mask']
    parameters=jnp.array(saved['parameters'])
    _,_,pa,pb=map(np.asarray,problem.critical(parameters))
    alpha=float(problem.alpha(parameters))
    coils=problem.finite_coils
    rng=np.random.default_rng(7654321)
    points=[rng.uniform([.5,-9],[14,9],(512,2)),rng.uniform([3.96,-4.55],[8.39,4.71],(256,2)),
            np.array([[.5,-9],[14,9],[.5,9],[14,-9],*np.asarray(problem.critical(parameters)[:2])])]
    # Coil interiors, edges, vertices, and both sides of the near/far switch.
    for r0,r1,z0,z1 in coils.rectangles:
        rc,zc=(r0+r1)/2,(z0+z1)/2
        gap=.5*max(r1-r0,z1-z0)
        points.append(np.array([[rc,zc],[r0,z0],[r1,z1],[r0,zc],[r1,zc],
                                [r1+gap-1e-7,zc],[r1+gap+1e-7,zc]]))
    # Probe both sides of saved support-cell boundaries within the chamber.
    mask=problem.mask
    k,i=np.where(mask[:,1:]!=mask[:,:-1])
    if len(i):
        select=np.linspace(0,len(i)-1,min(100,len(i)),dtype=int)
        for eps in [-1e-8,1e-8]:
            points.append(np.column_stack((3+(i[select]+1)*problem.dr+eps,
                                          -5.3+(k[select]+.5)*problem.dz)))
    points=np.concatenate(points)
    evaluate_self=jax.jit(jax.vmap(problem.self_flux,in_axes=(None,0)))
    selfpsi=np.asarray(evaluate_self(parameters,jnp.array(points)))*FLUX
    psi=selfpsi+coils.basis(points)@(np.asarray(problem.coils(parameters))*ICOIL)
    active=problem.mask_at(points)&(psi>float(pb)*FLUX)
    _,primitive=taper((float(pa)-psi/FLUX)/(float(pa)-float(pb)),.15,np)
    g=np.where(active,FB-alpha*FLUX*(float(pa)-float(pb))*primitive,FB)
    np.savetxt(out/'points.txt',points,fmt='%.17g')
    subprocess.run([str(args.executable.resolve()),'points',str(out/'points.txt'),str(out/'cpp_values.txt')],check=True)
    values=np.loadtxt(out/'cpp_values.txt')
    psi_error=float(np.max(abs(psi-values[:,2])))
    g_error=float(np.max(abs(g-values[:,3])))
    assert psi_error<2e-8,psi_error
    assert g_error<2e-9,g_error
    assert np.max(abs(values[:,4]+values[:,2]))==0
    # Compare cached plotted fields too, tying the exported evaluator to the figure.
    fields=np.load(snapshot/'fields.npz',allow_pickle=False)
    zi=rng.integers(0,len(fields['z']),128);ri=rng.integers(0,len(fields['r']),128)
    cache_points=np.column_stack((fields['r'][ri],fields['z'][zi]))
    np.savetxt(out/'cache_points.txt',cache_points,fmt='%.17g')
    subprocess.run([str(args.executable.resolve()),'points',str(out/'cache_points.txt'),str(out/'cache_values.txt')],check=True)
    cache=np.loadtxt(out/'cache_values.txt')
    cached_psi_error=float(np.max(abs(cache[:,2]-fields['psi_Vs'][zi,ri])))
    cached_g_error=float(np.max(abs(cache[:,3]-fields['B_R_Z_phi_T'][zi,ri,2]*cache_points[:,0])))
    assert cached_psi_error<2e-8,cached_psi_error
    assert cached_g_error<2e-9,cached_g_error
    # Full production-size grid and all toroidal planes, including the periodic face.
    subprocess.run([str(args.executable.resolve()),'grid',str(out/'grid100x02x200'),'100','2','200'],check=True)
    def read_vector(path):
        with path.open() as file:
            count=int(file.readline());data=np.loadtxt(file)
        assert len(data)==count
        return data
    psi_grid=read_vector(out/'grid100x02x200/vecpsi_grid100x02x200.txt').reshape(201,2,101)
    g_grid=read_vector(out/'grid100x02x200/vecg_grid100x02x200.txt').reshape(200,3,100)
    assert np.array_equal(psi_grid[:,0],psi_grid[:,1])
    assert np.array_equal(g_grid[:,0],g_grid[:,1]) and np.array_equal(g_grid[:,0],g_grid[:,2])
    # Check independently chosen grid nodes and centers against Python, avoiding
    # self-consistency-only validation of the C++ grid implementation.
    k=rng.integers(0,201,32);i=rng.integers(0,101,32)
    q=np.column_stack((3.05+i*(9.95-3.05)/100,-5.95+k*11.9/200))
    pypsi=np.asarray(evaluate_self(parameters,jnp.array(q)))*FLUX+coils.basis(q)@(np.asarray(problem.coils(parameters))*ICOIL)
    grid_error=float(np.max(abs(psi_grid[k,0,i]+pypsi)))
    assert grid_error<2e-8,grid_error
    k=rng.integers(0,200,32);i=rng.integers(0,100,32)
    q=np.column_stack((3.05+(i+.5)*(9.95-3.05)/100,-5.95+(k+.5)*11.9/200))
    pypsi=np.asarray(evaluate_self(parameters,jnp.array(q)))*FLUX+coils.basis(q)@(np.asarray(problem.coils(parameters))*ICOIL)
    _,primitive=taper((float(pa)-pypsi/FLUX)/(float(pa)-float(pb)),.15,np)
    pyg=np.where(problem.mask_at(q)&(pypsi>float(pb)*FLUX),FB-alpha*FLUX*(float(pa)-float(pb))*primitive,FB)
    grid_g_error=float(np.max(abs(g_grid[k,0,i]-pyg)))
    assert grid_g_error<2e-9,grid_g_error
    report=dict(passed=True,points=len(points),psi_max_abs_error_Vs=psi_error,g_max_abs_error_Tm=g_error,
        cached_plot_psi_max_abs_error_Vs=cached_psi_error,cached_plot_g_max_abs_error_Tm=cached_g_error,
        grid_psi_max_abs_error_Vs=grid_error,grid_g_max_abs_error_Tm=grid_g_error,
        source_axis_flux_Vs=float(pa)*FLUX,source_xpoint_flux_Vs=float(pb)*FLUX,
        all_toroidal_planes_identical=True,grid_counts=[psi_grid.size,g_grid.size],
        seconds=time.perf_counter()-started,checkpoint_sha256=hashlib.sha256((snapshot/'checkpoint.npz').read_bytes()).hexdigest(),
        executable_sha256=hashlib.sha256(args.executable.read_bytes()).hexdigest(),jax_device=str(jax.devices()[0]))
    (out/'validation.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2),flush=True)

if __name__=='__main__':main()
