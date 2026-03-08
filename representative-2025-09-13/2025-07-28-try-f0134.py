load("2025-07-23-representative.py")
xm_nox()
gb_verbose(1)

ring=PolynomialRing(QQ,'x,dx,y,dy,a,b')
x,dx,y,dy,a,b=ring.gens()
dic={'x':x,'dx':dx,'y':y,'dy':dy,'a':a,'b':b}
#cont=f0134_contiguity([a,b],[a-1,b],a.parent())  # for test
f=representative([[1,0],[0,1]],[0,0],f0134_contiguity,vars=[a,b])
save(f,"tmp-f0134-representative.sb")

