import itertools
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

# ---------- 1. Diamond unit cell ----------
corners = np.array(list(itertools.product([0, 1], repeat=3)), float)
faces = np.array([[.5, .5, 0], [.5, .5, 1], [.5, 0, .5], [.5, 1, .5], [0, .5, .5], [1, .5, .5]])
inner = np.array([[.25, .25, .25], [.75, .75, .25], [.75, .25, .75], [.25, .75, .75]])
atoms = np.vstack([corners, faces, inner])
bonds = [(a, b) for a, b in itertools.combinations(atoms, 2)
         if abs(np.linalg.norm(a - b) - np.sqrt(3) / 4) < 1e-6]

def draw_cell(ax):
    for a, b in itertools.combinations(corners, 2):          # cube edges
        if np.isclose(np.linalg.norm(a - b), 1):
            ax.plot(*zip(a, b), c='gray', lw=.6)
    for a, b in bonds:                                       # C-C bonds
        ax.plot(*zip(a, b), c='k', lw=1)
    ax.scatter(*np.vstack([corners, faces]).T, c='lightgray', edgecolors='k', s=50)
    ax.scatter(*inner.T, c='tab:orange', edgecolors='k', s=60)
    ax.set_box_aspect((1, 1, 1))

# ---------- 2. NV axes: the four <111> directions ----------
axes111 = np.array([[1, 1, 1], [-1, -1, 1], [-1, 1, -1], [1, -1, -1]]) / np.sqrt(3)

def to_lab(R, vecs):
    """R rows = lab x,y,z in crystal coords. Returns the n_z>0 member; -n is drawn too."""
    out = vecs @ R.T
    return out * np.sign(out[:, 2])[:, None]

orient = {
    '(100)': dict(
        R=np.array([[1, 1, 0] / np.sqrt(2), [1, -1, 0] / np.sqrt(2), [0, 0, 1]]),
        plane=[[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]]),
    '(111)': dict(
        R=np.array([[-1, 1, 0] / np.sqrt(2), [-1, -1, 2] / np.sqrt(6), [1, 1, 1] / np.sqrt(3)]),
        plane=[[1, 0, 0], [0, 1, 0], [0, 0, 1]]),
}

cols = ['tab:red', 'tab:blue', 'tab:green', 'tab:purple']

# ---------- 3. Figure 1: cell + both planes ----------
fig = plt.figure(figsize=(11, 5))
for i, (name, d) in enumerate(orient.items(), 1):
    ax = fig.add_subplot(1, 2, i, projection='3d')
    draw_cell(ax)
    ax.add_collection3d(Poly3DCollection([d['plane']], alpha=.35, fc='cyan', ec='b'))
    ax.set_title(f'Diamond cubic cell with {name} plane')
    ax.set_xlabel('[100]'); ax.set_ylabel('[010]'); ax.set_zlabel('[001]')
plt.tight_layout()
plt.savefig('diamond_planes.png', dpi=150)

# ---------- 4. Figures 2 & 3: one per orientation ----------
for name, d in orient.items():
    R = d['R']
    lab = to_lab(R, axes111)
    fig = plt.figure(figsize=(12, 5.5))

    # (a) crystal frame: cell, plane, NV axes drawn from cell centre
    ax = fig.add_subplot(1, 2, 1, projection='3d')
    draw_cell(ax)
    ax.add_collection3d(Poly3DCollection([d['plane']], alpha=.3, fc='cyan', ec='b'))
    c = np.array([.5, .5, .5])
    for v, col, k in zip(axes111, cols, range(1, 5)):
        ax.quiver(*c, *(.5 * v), color=col, lw=2, arrow_length_ratio=.15)
        ax.quiver(*c, *(-.5 * v), color=col, lw=2, ls='--', arrow_length_ratio=.15)
    ax.set_title(f'{name}: crystal frame, NV axes = <111>')
    ax.set_xlabel('[100]'); ax.set_ylabel('[010]'); ax.set_zlabel('[001]')

    # (b) lab frame: cut plate in x-y plane, z = surface normal
    ax = fig.add_subplot(1, 2, 2, projection='3d')
    L, t = 1.2, 0.08
    v = np.array(list(itertools.product([-L, L], [-L, L], [-t, t])))
    idx = [[0,1,3,2],[4,5,7,6],[0,1,5,4],[2,3,7,6],[0,2,6,4],[1,3,7,5]]
    ax.add_collection3d(Poly3DCollection([v[i] for i in idx], alpha=.12, fc='cyan', ec='gray'))
    for k, (n, col) in enumerate(zip(lab, cols), 1):
        ax.quiver(0, 0, 0, *n, color=col, lw=2.5, arrow_length_ratio=.12,
                  label=f'+NV{k}: ({n[0]:+.3f}, {n[1]:+.3f}, {n[2]:+.3f})')
        ax.quiver(0, 0, 0, *(-n), color=col, lw=2.5, ls='--', arrow_length_ratio=.12,
                  label=f'-NV{k}: ({-n[0]:+.3f}, {-n[1]:+.3f}, {-n[2]:+.3f})')
        for sgn in (1, -1):                                          # in-plane projections
            ax.plot([0, sgn * n[0]], [0, sgn * n[1]], [0, 0], c=col, ls=':')
    ax.quiver(0, 0, -1.4, 0, 0, 2.8, color='k', lw=1, arrow_length_ratio=.05)
    ax.text(0, 0, 1.5, 'surface normal z')
    ax.set_xlim(-L, L); ax.set_ylim(-L, L); ax.set_zlim(-1.5, 1.5)
    ax.set_box_aspect((1, 1, 1))
    ax.set_xlabel('x'); ax.set_ylabel('y'); ax.set_zlabel('z')
    ax.set_title(f'{name} diamond: lab frame (cut plane = x-y)')
    ax.legend(loc='upper left', fontsize=6, ncol=2)

    ang = np.degrees(np.arccos(np.clip(lab[:, 2], -1, 1)))
    print(name, 'angles to normal (deg):', np.round(ang, 2))
    print(np.round(np.vstack([lab, -lab]), 4))
    plt.tight_layout()
    plt.savefig(f'nv_{name.strip("()")}.png', dpi=150)

plt.show()
