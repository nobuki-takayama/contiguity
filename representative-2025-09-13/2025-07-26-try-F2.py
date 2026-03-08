load("2025-07-23-representative.py")
ring=PolynomialRing(QQ,'x,dx,y,dy,a,b,bp,c,cp')
x,dx,y,dy,a,b,bp,c,cp=ring.gens()
dic={'x':x,'dx':dx,'y':y,'dy':dy,'a':a,'b':b,'bp':bp,'c':c,'cp':cp}
f=representative([[1,0,0,0,0],[0,1,0,0,0],[0,0,1,0,0],[0,0,0,1,0],[0,0,0,0,1]],[0,0,0,0,0],f2_contiguity,vars=[a,b,bp,c,cp])
save(f,"tmp-F2-representative.sb")

