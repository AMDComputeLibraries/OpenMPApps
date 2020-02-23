#!/usr/bin/env python2
"""

make_input.py

Makes SIZE and data.rea files with the specified polynomial order and number of elements

Usage:
    $ ./make_input.py --lx1 <polynomial order> --lelt <number of elements>


"""

def nb_make_input(lx1, lelt):

    size_template = \
"""C     Dimension file to be included

      parameter (ldim=3)                      ! dimension
      parameter (lx1={lx1},ly1=lx1,lz1=lx1)      ! polynomial order

      parameter (lp = 10)                     ! max number of processors
      parameter (lelt={lelt})                    ! max number of elements, per proc

      parameter (lelg=lelt*lp)                ! max total elements in a test

      parameter (ldimt=1,ldimt1=ldimt+1)      ! used in 'include' files


      common /dimn/ nelx,nely,nelz,nelt       ! local element common block
     $            , nx1,ny1,nz1,ndim,nfield,nid
"""
    rea_template = \
""".true.    = ifbrick                  ! brick or linear geometry
{lelt} {lelt} 1 = iel0,ielN(per proc),stride ! range of number of elements per proc.
{lx1} {lx1}  1 = nx0,nxN,stride           ! poly. order range for nx1
0    0  0 = npx,npy,npz              ! np distrb, if != np, nekbone handle
0    0  0 = mx,my,mz                 ! nelt distrb, if != nelt, nekbone handle
"""

    with open('SIZE', 'w') as f:
        f.write(size_template.format(lx1=lx1, lelt=lelt))

    with open('data.rea', 'w') as f:
        f.write(rea_template.format(lx1=lx1, lelt=lelt))


if __name__ == '__main__':

    import argparse
    p = argparse.ArgumentParser(description='Create SIZE and data.rea files with the given parameters')
    p.add_argument('--lx1', type=int, default=16, help='the polynomial order')
    p.add_argument('--lelt', type=int, default=256, help='the number of elements')
    args = p.parse_args()

    nb_make_input(lx1=args.lx1, lelt=args.lelt)
