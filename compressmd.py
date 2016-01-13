#!/usr/bin/env python
import MDAnalysis, argparse, math, tempfile, shutil, os, sys
import numpy as np
from PIL import Image


parser = argparse.ArgumentParser(description='Compress a MD trajectory')
parser.add_argument('topo',metavar='topology_file')
parser.add_argument('traj',metavar='trajectory_file')
parser.add_argument('--out',metavar='output_file',default="output.m4v")
parser.add_argument("--selection",default="protein",required=False)

args = parser.parse_args()

model = MDAnalysis.Universe(args.topo,args.traj)
selection = model.select_atoms(args.selection)
n = selection.n_atoms

#make square-ish images -- this can be optimized to reduce dead pixels
#these also have to be even..
h = int(math.ceil(math.sqrt(n)))
w = h

#figure out overall range for x,y,z
maxvals = np.array([float("-inf")]*3)
minvals = np.array([float("inf")]*3)

for ts in model.trajectory:
	for c in selection.coordinates():
		for i in xrange(3):
			maxvals[i] = max(maxvals[i],c[i])
			minvals[i] = min(minvals[i],c[i])

ranges = maxvals-minvals;

print n,model.trajectory.n_frames
print minvals, ranges
#create an image for each frame in a tempdir
tmpdir = tempfile.mkdtemp()
for ts in model.trajectory:
	img = Image.new("RGB",(w,h))
	pix = img.load()
	#todo: support more than 8bits per coordinate
	for (ci,c) in enumerate(selection.coordinates()):
		i = ci/h
		j = ci%h
		pt = [0,0,0]
		for dim in xrange(3):
			pt[dim] = int(round(255*(c[dim]-minvals[dim])/ranges[dim]))
		pix[i,j] = tuple(pt)
	img.save('%s/%d.png' % (tmpdir, ts.frame))

os.system('ffmpeg  -i %s/%%d.png -c:v libx264  -preset veryslow -qp 0 %s' % (tmpdir, args.out))
shutil.rmtree(tmpdir)
