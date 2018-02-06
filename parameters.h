      COMMON/LOCAL/TILTL(nrs),ROLLL(nrs),TWISTL(nrs),SLY(nrs),SLX(nrs),
     1     DZL(nrs),LNDEX,DIETLT(nrs),DIEROL(nrs),PROTW(nrs),BUCAN(nrs),
     2     SLZ(nrs),TWISTH(nrs),OPENAN(nrs),OPENDS(nrs),ANGBN9(nrs),
     3     OPENC1(nrs),ANGBN1(nrs),APROP(nrs),DPROP(nrs),DZBUCK(nrs),
     4     DOPEN(nrs),CCDT(4000,4000),PPDT(4000,4000),OVRLP(nrs)
      common/locals/tilts(nrs),rolls(nrs),twists(nrs),slys(nrs),
     1     slxs(nrs),dzls(nrs)
      COMMON /DISPL/DXLOC(nrs),DYLOC(nrs),DZLOC(nrs),KKTI,
     1      BSCENTR(nrs,3),HORIG(nrs,3)
      COMMON/GLOB/TILTG(nrs),ROLLG(nrs),TWISTG(nrs),DXG(nrs),DYG(nrs),
     1           HG(nrs)
        common /bpwedg/bptilt(nrs),bproll(nrs),bptwst(nrs),bpshft(nrs),
     1                 bpslid(nrs),bprise(nrs),npradd,npadd1(nrs),
     2                 npadd2(nrs),badd1(nrs),badd2(nrs)
      COMMON /ZSTRAXS/ZAXISH(3,nrs),gaxish(3)
	double precision zaxish,gaxish
	character*1 badd1,badd2

