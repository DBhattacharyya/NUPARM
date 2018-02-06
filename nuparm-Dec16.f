
!==============================================================================
!  CORRECTED ON 11.12.2013  (Sankar, Manas, Pavan, Debasish)
!  ALL SUBROUTINES (SURFACE & OVERLAP & AXIS) CORRECTED AND UPDATED
!==============================================================================
C
C Please copy the AdeVariants.name, GuaVariants.name, CytVariants.name,
C UraVariants.name and surface.xyz files into a suitable directory, modify 
C lines containing open statement to open these files and recompile before
C using this program.  In this present form the program expects these files
C in /usr/local/bin, hence, if you can move the above files there you may
C not need to modify the program or recompile.
C



      common/radii/surf(50000,3),rad(50),istp(50),ienp(50),nrad
       COMMON /REPLY/ANSOR,ANSC1,ANSBPN,ANSDB,ANSrnt,ANSHO,ansnrm,answw,
     1    anspp,anstor,ansovl,anshlx

        dimension vv(10),xbp(100),ybp(100),zbp(100),el1(100),el2(100)
     1       ,el3(100)
        character *132 line,xis
        character*80  arguments,foutname
        character*40 valp, val1(20),input,filenm,fpdbfl,fprmfl,hlxfl
        character*40 param
        character*1 ansdb,ansbpn,ansnrm,ansor,ansc1,ansrnt,ansho,answw,
     1 anssng,anszp,anstrj,anspp,anstor,angl,ansovl,anshlx,axisb
        character*4 anorin
        logical unknopt
        character*8 cvar(999),cx
        real variable(99)

        ansdb='Y'
        ansbpn='Y'
        ansnrm='Y'
        anssng='N'
        ansor='N'
        ansc1='N'
        ansrnt='N'
        answw='N'
        anstor='N'
        input='        '
        anorin='    '
        anspp='N'
        anszp='N'
        ansovl='N'
        anshlx='N'
        hlxfl='        '
        anstrj='N'
        angl='N'
        lc=0
        axisb='N'

        npass=0
	nargument=iargc()
	if(nargument.eq.0) call errstop()
        do i=1,nargument
           call getarg(i,val1(i))
        enddo
        i=1
	do while(i.le.nargument)
	  if(val1(i)(1:1).eq.'-') then
            valp=val1(i)(2:20)
            unknopt=.true.
            if(valp(1:6).eq.'notpdb') then
              write(6,*) 'Please run NUPARM in interactive mode'
              stop
            elseif(valp(1:9).eq.'notdouble') then
              ansdb = 'N'
              unknopt=.false.
            elseif(valp(1:3).eq.'trj'.or.valp(1:4).eq.'traj'.or.
     1              valp(1:3).eq.'dcd'.or.valp(1:3).eq.'DCD') then
              anstrj='Y'
              unknopt=.false.
            elseif(valp(1:5).eq.'bpinf') then
              ansbpn = 'N'
              input = val1(i+1)
              i = i+1
              unknopt=.false.
            elseif(valp(1:3).eq.'hlx') then
              anshlx = 'Y'
              hlxfl = val1(i+1)
              i = i+1
              unknopt=.false.
            elseif(valp(1:5).eq.'lsfit') then
              ansnrm = 'N'
              n=rand(1)
              unknopt=.false.
            elseif(valp(1:2).eq.'cg') then
              ansor = 'Y'
              unknopt=.false.
            elseif(valp(1:4).eq.'c1c1') then
              ansc1 = 'Y'
              unknopt=.false.
              unknopt=.false.
            elseif(valp(1:6).eq.'single') then
              anssng = 'Y'
              unknopt=.false.
            elseif(valp(1:6).eq.'orient'.or.valp(1:8).eq.'reorient')then
              ansrnt = 'Y'
              anorin = val1(i+1)
              i=i+1
              unknopt=.false.
            elseif(valp(1:2).eq.'pp') then
              anspp = 'Y'
              unknopt=.false.
            elseif(valp(1:7).eq.'torsion'.or.valp(1:3).eq.'tor') then
              anstor = 'Y'
              unknopt=.false.
            elseif(valp(1:8).eq.'overlap'.or.valp(1:3).eq.'ovl') then
              ansovl = 'Y'
              call readsurf
              unknopt=.false.
            elseif(valp(1:2).eq.'ww'.or.valp(1:2).eq.'WW') then
              answw = 'Y'
              unknopt = .false.
            elseif(valp(1:3).eq.'par'.or.valp(1:3).eq.'PAR') then
              unknopt=.false.
              param=val1(i+1)
              i=i+1
            elseif(valp(1:3).eq.'out'.or.valp(1:3).eq.'OUT') then
              unknopt=.false.
              foutname=val1(i+1)
              i=i+1
            elseif(valp(1:4).eq.'AXIS' .or. valp(1:4).eq.'axis') then
              axisb = 'Y'
              unknopt=.false.
            endif
              if(unknopt) call errstop()
            else
            filenm = val1(i)
          endif
          i=i+1
        enddo   
        if(param.eq.'pp') anspp='Y'
        if(param.eq.'ovl') then
           ansovl='Y'
           call readsurf
        endif
        if(param.eq.'zp'.or.param.eq.'Zp'.or.param.eq.'ZP') anszp='Y'
        if(anszp.eq.'Y') anspp = 'Y'
        if(param.eq.'ang') angl='Y'

C        if(anstrj.eq.'N'.and.anorin.ne.'    ') then
        if(anstrj.eq.'N') then
         if(anorin.ne.'    ') write(*,*) 
     1                  'Molecule to be oriented using',anorin
          call nup(filenm,npass,input,anorin,hlxfl)
          if(axisb .eq. 'Y') then 
!           write (*,*) 'the value is :',axisb,' PDB :',filenm
           call axisal(filenm)
          endif
!          rewind(unit=11)
!          do while(1 .eq. 1)
!          read(11,88,end=100) xis
!88        format(a125)
!          enddo
!100       continue
          stop
        else
 
	  write(*,*) 'Reading trajector file: ',filenm
	  write(*,*) 'Writing parameters in: ',foutname
	  write(*,*) 'Calculating ',param,' for the snapshots'
          open(unit=17,file=foutname)
          if(param(1:5).eq.'alpha'.or.param(1:4).eq.'beta'.or.param(1:4)
     1 .eq.'gamma'.or.param(1:5).eq.'delta'.or.param(1:3).eq.'eps'.or.
     2 param(1:4).eq.'zeta'.or.param(1:3).eq.'chi'.or.param(1:3).eq.
     3 'amp'.or.param(1:5).eq.'phase'.or.param(1:3).eq.'eta'.or.
     4 param(1:5).eq.'theta') npass=1
          if(param(1:4).eq.'roll'.or.param(1:4).eq.'tilt'.or.param(1:5).
     1 eq.'twist'.or.param(1:5).eq.'shift'.or.param(1:5).eq.'slide'.or.
     2 param(1:4).eq.'rise'.or.param(1:3).eq.'cup') lc=1
          open(unit=12,file=filenm)
c	  call getlog(line)
          idpid=getpid()
          if(idpid.le.9) then
            write(fpdbfl,101) idpid
            write(fprmfl,1001) idpid
          elseif(idpid.le.99) then
            write(fpdbfl,102) idpid
            write(fprmfl,1002) idpid
          elseif(idpid.le.999) then
            write(fpdbfl,103) idpid
            write(fprmfl,1003) idpid
          elseif(idpid.le.9999) then
            write(fpdbfl,104) idpid
            write(fprmfl,1004) idpid
          elseif(idpid.le.99999) then
            write(fpdbfl,105) idpid
            write(fprmfl,1005) idpid
          else
             write(fpdbfl,106) idpid
            write(fprmfl,1006) idpid
          endif
101       format('/tmp/xx',i1,'.pdb')
102       format('/tmp/xx',i2,'.pdb')
103       format('/tmp/xx',i3,'.pdb')
104       format('/tmp/xx',i4,'.pdb')
105       format('/tmp/xx',i5,'.pdb')
106       format('/tmp/xx',i6,'.pdb')
1001       format('/tmp/xx',i1,'.prm')
1002       format('/tmp/xx',i2,'.prm')
1003       format('/tmp/xx',i3,'.prm')
1004       format('/tmp/xx',i4,'.prm')
1005       format('/tmp/xx',i5,'.prm')
1006       format('/tmp/xx',i6,'.prm')

c	  i=index(line,' ')
c	  filecr='/tmp/'//line(1:i-1)//'.pdb'
c	  filepr='/tmp/'//line(1:i-1)//'.prm'
	  open(unit=3,file=fpdbfl)
          i=0
          j=0
          kount=1
          do while(i.eq.0)
            read(12,1,end=99) line
            if(line(1:3).ne.'END') then
              write(3,1) line
            else
              close(unit=3)
              call nup(fpdbfl,npass,input,anorin,hlxfl)
              open(unit=10,file=fprmfl)
              kounter=0
              do while(j.le.0) 
                read(10,1,end=98) line
C
C   Reading Local Doublet Parameters
C
                if(param(1:4).eq.'tilt') then
                  if(line(1:2).eq.'LC') then
                    kounter=kounter+1
                    read(line,2) cvar(kounter)
                  endif
                elseif (param(1:4).eq.'roll') then
                  if(line(1:2).eq.'LC') then
                    kounter=kounter+1
                   read(line,2) cx,cvar(kounter)
                  endif
                elseif(param(1:5).eq.'twist') then
                  if(line(1:2).eq.'LC') then
                    kounter=kounter+1
                    read(line,2) cx,cx,cvar(kounter)
                  endif
                elseif(param(1:5).eq.'shift') then
                  if(line(1:2).eq.'LC') then
                    kounter=kounter+1
                    read(line,2)cx,cx,cx,cvar(kounter)
                  endif
                elseif(param(1:5).eq.'slide') then
                  if(line(1:2).eq.'LC') then
                    kounter=kounter+1
                    read(line,2)cx,cx,cx,cx,cvar(kounter)
                  endif
                elseif(param(1:4).eq.'rise') then
                  if(line(1:2).eq.'LC') then
                    kounter=kounter+1
                    read(line,2)cx,cx,cx,cx,cx,cvar(kounter)
                  endif
                elseif(param(1:3).eq.'cup') then
                  if(line(1:2).eq.'LC') then
                    kounter=kounter+1
                    read(line,2)cx,cx,cx,cx,cx,cx,cvar(kounter)
                  endif
		elseif(param(1:5).eq.'ovl') then
		  if(line(1:2).eq.'LC') then
		    kounter=kounter+1
		    read(line,2)cx,cx,cx,cx,cx,cx,cx,cvar(kounter)
		  endif
C 
C Parsing Base Pair Parameters
C
                elseif(param(1:6).eq.'buckle') then
                    if(line(1:2).eq.'BL') then
                      kounter=kounter+1
                      read(line,2)cvar(kounter)
                    endif
                elseif(param(1:4).eq.'open') then
                   if(line(1:2).eq.'BL') then
                     kounter=kounter+1
                     read(line,2)cx,cvar(kounter)
                   endif
                elseif(param(1:4).eq.'prop') then
                   if(line(1:2).eq.'BL') then
                     kounter=kounter+1
                     read(line,2)cx,cx,cvar(kounter)
                   endif
                elseif(param(1:4).eq.'stag') then
                  if(line(1:2).eq.'BL') then
                     kounter=kounter+1
                     read(line,2)cx,cx,cx,cvar(kounter)
                  endif
                elseif(param(1:5).eq.'shear') then
                   if(line(1:2).eq.'BL')then
                     kounter=kounter+1
                     read(line,2)cx,cx,cx,cx,cvar(kounter)
                   endif
                elseif(param(1:7).eq.'stretch') then
                   if(line(1:2).eq.'BL') then
                     kounter=kounter+1
                     read(line,2)cx,cx,cx,cx,cx,cvar(kounter)
                   endif
C
C   Reading Torsion Angles
C
                elseif(param(1:5).eq.'alpha') then
                    if(line(1:2).eq.'TR') then
                       kounter=kounter+1
                       read(line,4) variable(kounter)
                    endif
                elseif(param(1:4).eq.'beta') then
                  if(line(1:2).eq.'TR') then
                    kounter=kounter+1
                    read(line,4) x,variable(kounter)
                  endif
                elseif(param(1:4).eq.'gamma') then
                  if(line(1:2).eq.'TR') then
                    kounter=kounter+1
                    read(line,4) x,x,variable(kounter)
                  endif
                elseif(param(1:5).eq.'delta') then
                  if(line(1:2).eq.'TR') then
                      kounter=kounter+1
                     read(line,4) x,x,x,variable(kounter)
                  endif
                elseif(param(1:3).eq.'eps') then
                  if(line(1:2).eq.'TR') then
                       kounter=kounter+1
                       read(line,4) x,x,x,x,variable(kounter)
                  endif
                elseif(param(1:4).eq.'zeta') then
                  if(line(1:2).eq.'TR') then
                       kounter=kounter+1
                       read(line,4) x,x,x,x,x,variable(kounter)
                  endif
                elseif(param(1:3).eq.'chi') then
                  if(line(1:2).eq.'TR') then
                       kounter=kounter+1
                       read(line,4) x,x,x,x,x,x,variable(kounter)
                  endif
                elseif(param(1:3).eq.'eta') then
                  if(line(1:2).eq.'TR') then
                       kounter=kounter+1
                       read(line,4) x,x,x,x,x,x,x,variable(kounter)
                  endif
                elseif(param(1:5).eq.'theta') then
                  if(line(1:2).eq.'TR') then
                       kounter=kounter+1
                       read(line,4) x,x,x,x,x,x,x,x,variable(kounter)
                  endif
C
C   Reading Sugar Puckher Phase and Amplitude
C
                elseif(param(1:3).eq.'amp') then
                  if(line(1:2).eq.'SG') then
                    kounter=kounter+1
                    read(line,5) variable(kounter)
                  endif
                elseif(param(1:5).eq.'phase') then
                  if(line(1:2).eq.'SG') then
                    kounter=kounter+1
                    read(line,5) x,variable(kounter)
                  endif
                elseif(param(1:2).eq.'zp'.or.param(1:2).eq.'ZP') then
                  if(line(4:13).eq.'Strand: 2') exit
                  if(line(1:2).eq.'ZP') then
                    kounter=kounter+1
                    read(line,7) variable(kounter)
                  endif
                endif
                if(param.eq.'pp'.and.anszp.ne.'Y') then
                  if(line(1:2).eq.'PP') then
                    write(17,6) kount,line(6:132)
                  endif
                endif
                if(angl.eq.'Y') then
                  if(line(8:16).eq.'Base/B.P.') then
                    read(10,1) line
                    read(10,1) line
                 read(line,9)xbp(1),ybp(1),zbp(1),el1(1),el2(1),el3(1)
c	write(*,*) xbp(1),ybp(1),zbp(1)
                    nr=1
                    do while(nr.lt.100)
                      read(10,1) line
                      nr=nr+1
                      if(line(4:13).eq.'END-TO-END') exit
                        read(line,9)xbp(nr),ybp(nr),zbp(nr),el1(nr),
     1                   el2(nr),el3(nr)
                                  
c                      endif
                    enddo
                  
            cost1= (el1(1)*el1(nr-1)+el2(1)*el2(nr-1)+el3(1)*el3(nr-1))/
     1            ((el1(1)*el1(1)+el2(1)*el2(1)+el3(1)*el3(1))**0.5 *
     2           (el1(nr-1)*el1(nr-1)+el2(nr-1)*el2(nr-1)+el3(nr-1)*
     3           el3(nr-1))**0.5)
                theta1=(180.00/3.14)*acos(cost1)
              dist1= sqrt((xbp(1)-xbp(nr-1))**2 +(ybp(1)-ybp(nr-1))**2 +
     1                    (zbp(1)-zbp(nr-1))**2)

                                   
            cost2=(el1(2)*el1(nr-2)+el2(2)*el2(nr-2)+el3(2)*el3(nr-2))/
     1              ((el1(2)*el1(2)+el2(2)*el2(2)+el3(2)*el3(2))**0.5 *
     2              (el1(nr-2)*el1(nr-2)+el2(nr-2)*el2(nr-2)+
     3              el3(nr-2)*el3(nr-2))**0.5)
                  theta2=(180.00/3.14)*acos(cost2)
                dist2= sqrt((xbp(2)-xbp(nr-2))**2 +(ybp(2)-ybp(nr-2))**2
     1                   +(zbp(2)-zbp(nr-2))**2)

            variable(1)=theta1
            variable(2)=theta2
            variable(3)=dist1
            variable(4)=dist2
c	write(*,*) 'Variables',variable(1),variable(2),variable(3),
c     1 variable(4)
            kounter=4 

                 endif
               endif
                    
              enddo
               
98            continue
	      close(unit=10)
              if(npass .ne. 1)then
                do k=1,kounter
                  variable(k)=0.0
                  if((cvar(k)(6:8).ne.'NAN').and.(cvar(k)(6:8).ne.'Nan')
     1.and.(cvar(k)(6:8).ne.'NaN').and.(cvar(k)(6:8).ne.'nan')) then
                  write(23,*)cvar(k),variable(k)
                    read(cvar(k),23) variable(k)
                  else
                    variable(k)=999.99
                  endif
                enddo
              endif
              if(lc.eq.1) then
                write(17,3) kount,(variable(k),k=1,kounter-1)
              else
                write(17,3) kount,(variable(k),k=1,kounter)
              endif
              open(unit=3,file=fpdbfl)
              kount=kount+1
            endif
          enddo
        endif
99      continue
        stop
1       format(a132)
2       format(18x,8a8)
3       format(i7,999f8.2)
4       format(10x,9f9.2)
5       format(54x,2f9.2)
6       format(i5,a)
7       format(90x,f9.2)
8       format('+ Processing frame no',i6)
9       format(3x,3f8.3,5x,3f10.4)
10      format(3f8.3)
23      format(f8.2)
        end

        subroutine nup(filecr,npass,input,anorin,hlxfl)
        
C 
C  This program compiles well with gnu/g77 or other standard Unix fortran 
C  compilers.  However, the "gfortran" compiler on RedHat Enterprise Linux
C  does not make correct executable due to its IARGC routine problem.
C
C
C   Adopted new method for calculation of Intra Base Pair parameters following 
C   S. Mukherjee, M. Bansal and D. Bhattacharyya (2006) J. Comp. Aided Mol. Des. 
C   vol. 20, 629-645.
C   It still calculates in the older way also.
C
C   Added additional functionality to calculate base normal in exact sense and
C   not following Least Squares Fit method, which is unstable at times. (2003).

C   Modified on March 15, 2003, to consider BASES in a different way
C   Modified on 21st. Aug, 1992, to accept AMBER-PDB file as such.
C
C  The NUPARM program has been developed by M. Bansal and D. Bhattacharyya 
C  of Molecular Biophysics Unit, Indian Institute of Science, Bangalore 
C  560012, India, so as to meet almost all the requirements stipulated at the 
C  EMBO Workshop on DNA Curvature and Bending, held at Cambridge, U.K., in
C  Sept. 1988.  The Nomenclature and Description of the DNA structural 
C  parameters follows the Cambridge Convention (as published in J. Mol. Biol.
C  vol. 205, 787-791 (1989), J. Biomol. Struct. Dynam, vol.6, 627-634 (1989)
C  EMBO J. vol. 8, 1-4 (1989)).  The definition of the local helix and wedge 
C  parameters are in terms of the local helix axis and the mean Z-axis
C  respectively for the doublet involved (as described in J. Biomol. Struct.
C  Dynam. vol. 6, 635-653 (1989)).  
C
C  If any bugs/errors are detected please inform the authors. mb@mbu.iisc.ernet.in 
C  or dhananjay.bhattacharyya@saha.ac.in
C
C-----------------------------------------------------------------------
	parameter ( natm = 999999 )
	parameter ( nrs = 9999 )
        include 'coordinates.h'
        include 'parameters.h'
      DIMENSION IBPAIR(nrs),IPAIR(4),IP(4),IN(4),DC(3),
     1           DCR(3),RMAT(3,3),TMP(3),TITLE1(20),CENTRS(nrs,3),
     2           ybase1(3),ybase2(3),xbase1(3),xbase2(3),zbase1(3),
     2           ybase3(3),xbase3(3),zbase3(3),zbase2(3),orig1(3),
     3           orig2(3),bm(3),xbs(2,3),ybs(2,3),param(6),orig3(3),
     4           orig4(3),DummyV(3),nbpadd(nrs,3),typeadd(nrs,3),
     5           bptadd(nrs,3),ybs1(2,3),ybs2(2,3)
      COMMON/ATOMS/IBASE1(10),IBASE2(10),IBASE3(10),IBASE4(10),NB,
     1ibapr(nrs)
      DOUBLE PRECISION XAXIS,YAXIS,ZAXIS,XM,YM,ZM1,ZM2,ZM3,Y1,Y2,Y3,
     1   Y12,Y22,Y32,X11,X12,X13,BMN,TL2,RL2,WTILT,WROLL,Y11,Y13,
     2   xbs,ybs,xbase1,ybase1,zbase1,xbase2,ybase2,zbase2,xbase3,
     3   ybase3,zbase3,DC,DCR1,DCR2,DCR3,ANGL,ZGX,ZGY,ZGZ,ybs1,ybs2
        COMMON /BASINF/BASE(nrs),pair(3,nrs),prtype(nrs)
      COMMON /GLOREO/DCRS1,DCRS2,DCRS3,ANGLE,DX1,DY1,INFBP(2,nrs)
C      COMMON /GLOREO/DCRS1,DCRS2,DCRS3,ANGLE,DX1,DY1,INFBP(2,nrs)
      COMMON /REPLY/ANSOR,ANSC1,ANSBPN,ANSDB,ANSrnt,ANSHO,ansnrm,answw,
     1   anspp,anstor,ansovl,anshlx
	logical stk,unknopt
	external matmul
       CHARACTER *80 FILENM,TITLE1,TITLE
       CHARACTER *40 FORMA,INPUT,valp,filecr,val1(20),filprt*20,hlxfl
       CHARACTER *4  ATOM,ANORIN,ATTYPE
       CHARACTER *3 BAS1,BAS2,BAS3,BAS4,BASE,type,typeadd,bptp
       CHARACTER *1 ANSrnt,ANSOR,ANSC1,ANSHO,ANSBPN,ANSPDB,ANSPP,ANSTOR,
     1ANSDB,ANSSNG,ANSCCYL,type1,type2,pairtype(2*nrs),pair,type3,type4,
     2 bptype(nrs),cystrns,ansnrm,bptadd,cystrns2,ansbp1,answw,ansovl,
     3 prtype,anshlx
C
       DATA MAXSEG,MAXAT/nrs,natm/
C --  This programme can handle a maximum of NRS residues and NATM atoms.
C
c        ANORIN = '    '
        CONV=180.0/3.141592654
        INOHP = 0
        NB = 1
        IND = 0
	nseladdl=0
C
c	write(*,*) 'Input coordinates from: ',filecr
       OPEN(UNIT=8,FILE='RUN.ANS')
       OPEN(UNIT=9,FILE='XYZAXES.OUT')
C
  3     FORMAT(/' All the questions & answers are written into a file
     1 "RUN.ANS" ')
c131      WRITE(6,4)
         WRITE(8,4)
  4     FORMAT(/' Type in the name of the INPUT DATA file:')
c	nargumnt = iargc()
c	if(nargumnt.eq.0) call errstop()
	anspdb = 'Y'
        if(npass.eq.1) anstor='Y'
c	do i=1,nargumnt
c	   call getarg(i,val1(i))
c	enddo
c	filenm = '/tmp/xx.pdb'
  5     FORMAT(A80)
         WRITE(8,'(1X,''COORDINATES TAKEN FROM :''/1X,A80)') filecr
         WRITE(9,'(1X,''COORDINATES TAKEN FROM :''/1X,A80)') filecr
        OPEN(UNIT=4,FILE=filecr,STATUS='OLD',ERR=131)
C         WRITE(6,3)
C
        CALL GETCOR(4,IERR,ANSPDB)
C
          IF(IERR.EQ.1) THEN
             STOP
          END IF
c	WRITE(6,81) 
 81	FORMAT(/'Do you want base normals by Cross-Products[Y]',$)
  89    FORMAT(/' Should the GLOBAL AXIS be fitted to:'/' 1. Local HELIX
     1 ORIGINS (Type O)'/' 2. Base-Pair CENTERS (Type C)'/' 3. Any Backb
     2one atom (Type the atom name)')
	WRITE(8,81) 
	WRITE(8,'(5X,A1)') ANSNRM
C-----------------------------------------------------------------------
C
        IF(ANSPDB.EQ.'Y') THEN
          REWIND 4
          IT = 0
          KKTI = 0
          DO 2000 KL=1,200
             READ(4,'(A80)') TITLE
	      IF(TITLE(1:6).EQ.'HEADER') THEN
	         KKTI = KKTI + 1
	         TITLE1(KKTI) = TITLE
             ELSE IF(TITLE(1:6).EQ.'COMPND') THEN
                 KKTI = KKTI + 1
                 TITLE1(KKTI) = TITLE
             ELSE IF(TITLE(1:6).EQ.'AUTHOR') THEN
                 KKTI = KKTI + 1
                 TITLE1(KKTI) = TITLE
             ELSE IF(TITLE(1:4).EQ.'JRNL') THEN
                 KKTI = KKTI + 1
                 TITLE1(KKTI) = TITLE
             END IF
             IF(TITLE(1:4).EQ.'ATOM') IT = 1
             IF(KKTI.GT.19) IT = 1
	      IF(IT.EQ.1) GO TO 2001
 2000	   CONTINUE
 2001	 CONTINUE
      END IF
      CLOSE(UNIT=4)
C
c      WRITE(6,91)
 91   FORMAT(/'Does the molecule have double stranded helix? [Y]',$)
      WRITE(8,91)
c      READ(5,'(A1)')ANSDB
      IF(ANSDB.NE.'n'.AND.ANSDB.NE.'N') THEN
	ANSDB = 'Y'
      ELSE
	ANSDB = 'N'
      END IF
      WRITE(8,'(5X,A1)')ANSDB
      NDB = NSEG/2
C
c      WRITE(6,92)
 92   FORMAT(/'Are the residues numbered according to the Brookhaven for
     1mat? [Y]',$)
      WRITE(8,92)
c      READ(5,'(A1)') ANSBPN
      IF(ANSBPN.NE.'n'.AND.ANSBPN.NE.'N') THEN
 	ANSBPN = 'Y'
      ELSE
	ANSBPN = 'N'
      END IF
      WRITE(8,'(5X,A1)') ANSBPN
      IF(ANSBPN.EQ.'Y'.AND.ANSDB.EQ.'Y') THEN
           KOUNT = -1
           DO 2002 KKK=1,NDB
              KOUNT = KOUNT + 2
              IBAPR(KOUNT) = KKK
              IBAPR(KOUNT + 1) = NSEG - KKK + 1
 2002	    CONTINUE
c	write(6,'(''a  '',2i6)') (ibapr(kkk),kkk=1,kount+1)
           NB = NSEG
	 if(nb.gt.3) nseladdl=1
       END IF
C
       IF(ANSBPN .NE. 'Y') THEN
C
c133        WRITE(6,93)
 93   FORMAT(/' Give the file-name containing base-pair numbering inform
     1ation (in 2I5 Format):')
          WRITE(8,93)
c          READ(5,'(A40)') INPUT
       WRITE(8,'(1X,''Base-pair residue numbers from:''/1X,A40)')INPUT
c       WRITE(6,'(1X,''Base-pair residue numbers from:''/1X,A40)')INPUT
          OPEN(UNIT=2,FILE=INPUT,STATUS='OLD',ERR=132)
          DO 2003  KB = 1,MAXSEG
             READ(2,10,END=1001) IB1,IB2,type,bptype(nb),nbpadd(kb,1),
     1 typeadd(kb,1),bptadd(kb,1),nbpadd(kb,2),typeadd(kb,2),
     2 bptadd(kb,2),nbpadd(kb,3),typeadd(kb,3),bptadd(kb,3)
  10       FORMAT(I5,1x,4(i5,1x,a3,1x,a1,3x))
  11 	format(2i5,4x,a1,' : ',a1)
             IBAPR(NB) = IB1
             IBAPR(NB + 1) = IB2
	     pairtype(nb)=type(1:1)
	     pairtype(nb+1)=type(3:3)
c	write(8,11)ibapr(nb),ibapr(nb+1),pairtype(nb),pairtype(nb+1)
c	write(6,11)ibapr(nb),nbpadd(kb,1),typeadd(kb,1)(1:1),typeadd(kb,1)(3:3)
c	write(6,11)ibapr(nb),nbpadd(kb,3),typeadd(kb,3)(1:1),typeadd(kb,3)(3:3)
	     if(answw.eq.'Y') then
C	       write(*,*) 'WWC Option selected',nb
	       pairtype(nb)='W'
	       pairtype(nb+1)='W'
	       bptype(nb)='C'
	       nbpadd(kb,1)=0
	       nbpadd(kb,2)=0
	       nbpadd(kb,3)=0
	     endif
             NB = NB + 2
 2003	   CONTINUE
 1001    CONTINUE
c	write(6,'(20(2x,a1,1x,a1))') (pairtype(kc),kc=1,kb)
	 if(nb.gt.3) then
	   nseladdl=1
	 else
	   ibapr3=ib1
	   ibapr4=ib2
	   bptp = type
	   cystrns2 = bptype(nb-2)
	   base(ibapr3)(1:1) = secnm(ibapr3)(1:1)
	   base(ibapr3)(2:2) = ':'
	   base(ibapr3)(3:3) = secnm(ibapr4)(1:1)
	   pair(1,ibapr3) = bptp(1:1)
	   pair(2,ibapr3) = bptp(3:3)
	   pair(3,ibapr3) = cystrns2
	 endif
c 	write(*,*) 'NB for this system',nb,ibapr3,ibapr4,bptp,cystrns2
c	write(*,*) base(ibapr3),'  ',cystrns2
c
c   General analysis of oligonucleotide has been requested.  It would pass
c   through all the modules
c
          NB = NB - 1
        END IF
c	write(6,'(''b '',2i6)') (ibapr(kkk),kkk=1,kount+1)
c
c  Check validity of base pairing information suitable for re-orientation
c
	if(ansrnt.eq.'Y'.or.ansrnt.eq.'y') then
	  do kbspr=4,nb,2
c	    write(*,*) ibapr(kbspr),ibapr(kbspr-3)
c	    write(*,*) ibapr(kbspr-1),ibapr(kbspr-2)
	    if(ibapr(kbspr).eq.0.or.(ibapr(kbspr).eq.ibapr(kbspr-3).and.
     1      ibapr(kbspr-1).eq.ibapr(kbspr-2))) then
C	      write(6,71) 
	      write(8,81)
	      ansrnt = 'N'
	      exit
	    endif
	  enddo
	endif
71	format(' Reorientation would fail in this molecule, attempting'
     1  ,' to calculate parameters without that option.')
C
        IF(ANSDB.NE.'Y'.AND.ANSBPN.EQ.'Y')THEN
            KOUNT = 1
           DO 2004 KK = 1,NSEG
              IBAPR(KOUNT) = KK
              IBAPR(KOUNT + 1) = 0
              KOUNT = KOUNT + 2
 2004	    CONTINUE
          NB = NSEG * 2
	  if(nb.gt.0) nseladdl=1
        END IF
c	write(6,'(''c '',2i6)') (ibapr(kkk),kkk=1,kount+1)
C
C-----------------------------------------------------------------------
C-- OPTION TO USE 'Centre of Gravity' AS ORIGIN IS SUPPRESSED.
C   CAN BE ACTIVATED, "AT YOUR PERIL", BY DE-COMMENTING NEXT SIX LINES.
        IF(ANSDB.EQ.'Y') THEN
c           WRITE(6,96)
 96   FORMAT(/' Should the base-pair ORIGIN be at the Centre of Gravity?
     1'/'(NOT RECOMMENDED!! Default is midpoint of C6--C8:) [N]',$)
           WRITE(8,96)
c           READ(5,'(A1)') ANSOR
	   IF(ANSOR.EQ.'y'.OR.ANSOR.EQ.'Y') THEN
	      ANSOR = 'Y'
 	   ELSE
	      ANSOR = 'N'
	   END IF
           WRITE(8,'(5X,A1)') ANSOR
        END IF
C-----------------------------------------------------------------------
C       IF THE MOLECULE IS SINGLE STRANDED THIS QUESTION IS OF NO USE
        IF (ANSDB .EQ. 'Y') THEN
c          WRITE(6,97)
 97   FORMAT(/' Should the Y-axis be along the line joining C1'' atoms?'
     1 /'(If "No",Y-axis is taken along the C6--C8 direction:) [N]',$)
c          READ(5,'(A1)') ANSC1
  	   IF(ANSC1.EQ.'y'.OR.ANSC1.EQ.'Y') THEN
	     ANSC1 = 'Y'
	   ELSE
	     ANSC1 = 'N'
	   END IF
          WRITE(8,97)
          WRITE(8,'(5X,A1)') ANSC1
        END IF
C
C 
C-----------------------------------------------------------------------
       NRUN =1
	IF(ANSDB.EQ.'Y') THEN
c         WRITE(6,'(/''Are Single Strand parameters required? [N]'',$)')
	  WRITE(8,'(/''Are Single Strand parameters required?(Y/N)'')')
c	  READ(5,'(A1)') ANSSNG
	  WRITE(8,'(5X,A1)') ANSSNG
	  IF(ANSSNG.EQ.'y'.OR.ANSSNG.EQ.'Y') NRUN=3
	END IF
C
       DO 700  NTIM = 1,NRUN
          NG = 1
       IF(NTIM .EQ. 1)NG = 3
C
       DO 2005 NGLOB = 1,NG
         LNDEX = 0
         NSTEPS = 0
	 stk=.false.
	if(nseladdl.gt.0) then
c
c  Calculation of Base Pair Step (Local Doublet) parameters begin 
C
         DO 600 IBP = 1,(NB-2),2
           NSTEPS = NSTEPS + 1
           IF(IBP.EQ.NB-1) INOHP = 1
           IBAPR1 = IBAPR(IBP)
           IBAPR2 = IBAPR(IBP + 1)
	   type1=pairtype(ibp)
	   type2=pairtype(ibp+1)
	   cystrns=bptype(ibp)
	    pair(1,nsteps) = pairtype(ibp)
	    pair(2,nsteps) = pairtype(ibp+1)
	    pair(3,nsteps) = bptype(ibp)
c	write(*,*) 'bp1',ibp,ibapr1,type1,' & bp 2',ibp+1,ibapr2,type2
c	write(8,13) ibapr1,ibapr2,ibp,type1,type2
13	format(' BasePair betn.',2i5,' PairNo.',i5,' Type',2a1)
           IF(NTIM .EQ. 2)IBAPR2 = 0
           IF(NTIM .EQ. 3) THEN
             IBAPR1 = IBAPR(NB - IBP + 1)
             IBAPR2 = 0
           END IF
C
           if(ibapr1.ne.0) then         ! bhatta March 31, 2012
           CALL PICATP(IBAPR1,IBAPR2,IPAIR,IBASE1,IBASE2,NB1,NB2,INOHP,
     1           BAS1,BAS2)
            BASE(NSTEPS)(1:1) = BAS1(1:1)
            BASE(NSTEPS)(2:2) = ':'
            BASE(NSTEPS)(3:3) = BAS2(1:1)
            INFBP(1,NSTEPS) = IBAPR1
            INFBP(2,NSTEPS) = IBAPR2
c	write(8,12) base(nsteps),ibapr1,ibapr2,nb1,nb2,pair(1,nsteps),
c     1  pair(2,nsteps),pairtype(ibapr1),pairtype(ibapr2)
12	format(' BasePair ',a3,' Betn.',2i5,' No.ofAtoms',2i5,4(1x,a1))
C
            IP(1) = IPAIR(2)
            IP(2) =  IPAIR(4)
            IN(1) = IPAIR(1)
            IN(2) = IPAIR(3)
C
            IBAPR3 = IBAPR(IBP+2)
            IBAPR4 = IBAPR(IBP+3)
c	write(*,*) 'IBP=',ibp
	    type3=pairtype(ibp+2)
	    type4=pairtype(ibp+3)
c	write(*,'(2x,i3,a1,1x,i3,a1)') (kc,pairtype(kc),kc=1,50)
c	write(*,*)'pairtype', ibapr4,type4,pairtype(ibapr4)
c	   type3=pairtype(ibp+2)
c	   type4=pairtype(ibp+3)
c	write(*,*) 'bp 3',ibp+2,ibapr3,type3, ' & bp 4',ibp+3,ibapr4,type4
	   cystrns2=bptype(ibp+2)
            IF(NTIM .EQ. 2)IBAPR4 = 0
            IF(NTIM .EQ. 3) THEN
              IBAPR3 = IBAPR(NB - IBP -1)
              IBAPR4 = 0
            END IF
C
           CALL PICATP(IBAPR3,IBAPR4,IPAIR,IBASE3,IBASE4,NB3,NB4,INOHP,
     1           BAS3,BAS4)
           BASE(NSTEPS+1)(1:1) = BAS3(1:1)
           BASE(NSTEPS+1)(2:2) = ':'
           BASE(NSTEPS+1)(3:3) = BAS4(1:1)
           INFBP(1,NSTEPS+1) = IBAPR3
           INFBP(2,NSTEPS+1) = IBAPR4
	   pair(1,nsteps+1) = pairtype(ibp+2)
	   pair(2,nsteps+1) = pairtype(ibp+3)
	   pair(3,nsteps+1) = bptype(ibp+2)
c	write(8,12) base(nsteps+1),ibapr3,ibapr4,nb3,nb4,
c     1       pair(1,nsteps+1),pair(2,nsteps+1)
c	write(8,12) base(nsteps+1),ibapr3,ibapr4,nb3,nb4,
c     1 pair(1,nsteps+1),pair(2,nsteps+1),pairtype(ibapr3),
c     2 pairtype(ibapr4)
           IP(3) = IPAIR(2)
           IP(4) = IPAIR(4)
           IN(3) = IPAIR(1)
           IN(4) = IPAIR(3)
C
C Calculation of Single Strand parameter between IBAPR1 and IBAPR3
C
	   if(ntim.gt.1) then
	     call findybs(ibapr1,xbase1,ybase1,zbase1,ang1,orig1,'W','C')
	     call findybs(ibapr3,xbase3,ybase3,zbase3,ang3,orig3,'W','C')
	     ang=0.0
C
C  Base Y-axis becomes IUPAC X-axis, Base X-axis becomes IUPAC Z-axis and 
C  Base Z-axis becomes IUPAC Y-axis 
C
            do kx=1,3
              ang=ang+ybase1(kx)*ybase3(kx)
              bm(kx)=orig3(kx)-orig1(kx)
              xbs(1,kx)=ybase1(kx)
              xbs(2,kx)=ybase3(kx)
              ybs(1,kx)=zbase1(kx)
c              ybs(2,kx)=-xbase3(kx)
              ybs(2,kx)=zbase3(kx)
            enddo
	write(9,*)'Xbs & Ybs1 b4 WEGDMbs',(xbs(1,Lc),ybs(1,Lc),Lc=1,3)
	write(9,*)'Xbs & Ybs2 b4 WEGDMbs',(xbs(2,Lc),ybs(2,Lc),Lc=1,3)
            call wegdmbs(xbs,ybs,bm,param,stk)
c	    write(6,*) 'S-S Parameters between',Ibapr1,ibapr3
c	    write(*,*) (param(kk),kk=1,6)
	    tilts(nsteps) = param(1)
	    rolls(nsteps) = param(2)
	    twists(nsteps) = param(3)
	    slxs(nsteps) = param(4)
	    slys(nsteps) = param(5)
	    dzls(nsteps) = param(6)
	  endif
C
C  Calculation of Intra Base Pair parameters following Wedge Definition
C
	if(ibapr2.ne.0) then
	    bptp = type1//':'//type2
c	write(*,*) nsteps,bptp
	    call getblprm(ibapr1,ibapr2,bptp,cystrns,
     2 bptilt(nsteps),bproll(nsteps),bptwst(nsteps),
     2 bpshft(nsteps),bpslid(nsteps),bprise(nsteps),stk,ybs1)    !Sukanya, 12 Nov, 2013
	endif
	if(ibapr4.ne.0) then
	  bptp = type3//':'//type4
c	write(*,*) nsteps,bptp
	    call getblprm(ibapr3,ibapr4,bptp,cystrns2,
     2 bptilt(nsteps+1),bproll(nsteps+1),bptwst(nsteps+1),
     2 bpshft(nsteps+1),bpslid(nsteps+1),bprise(nsteps+1),stk,ybs2)    !Sukanya, 12 Nov, 2013
	endif

	write(9,*) 'Type 1, 2, 3 & 4 ',ibapr1,ibapr2,ibapr3,ibapr4,bptp
C
C  Calculation of Inter Base Pair Parameters is done mostly in the 
C  Subroutine below
C
c	write(*,*)IBAPR1,IBAPR2,IBAPR3,IBAPR4,' FINDVEC' 
          CALL FINDVEC(IP,IN,NB1,NB2,NB3,NB4,NSTEPS,ANSOR,ANSC1,ANSrnt,
     1    NTIM,ibapr1,ibapr3,ybs1,ybs2,type1,type2,type3,type4,cystrns,
     2    cystrns2)    !Sukanya, 12 Nov, 2013
        else                          ! bhatta March 31, 2012
          lndex = lndex + 1
        endif
 600  CONTINUE
	endif
	if(ibapr4.ne.0) then
	    if(nseladdl.ne.0) bptp = type3//':'//type4
	    prtype(nsteps+1) = cystrns2
	    call getblprm(ibapr3,ibapr4,bptp,cystrns2,
     2 bptilt(nsteps+1),bproll(nsteps+1),bptwst(nsteps+1),
     2 bpshft(nsteps+1),bpslid(nsteps+1),bprise(nsteps+1),stk,ybs)
c 	  write(*,*) nsteps,ibapr3,ibapr4,bptp,' ',cystrns2,bptilt(nsteps+1)
	endif
		
c
c   Parameter for a single base pair has been requested.  Hence local doublet
c   parameter calculation is not done and only additional base pair parameter
c   calculation is taking place
c
	npradd=0
c	write(*,*) 'Additional Base Pairs',npradd,nsteps
	do 601 i=1,nsteps+1
	  do 601 k=1,3
	  if(nbpadd(i,k).ne.0) then
	    npradd=npradd+1
	    iadd=nsteps+1+npradd
	    pair(1,iadd) = typeadd(i,k)(1:1)
	    pair(2,iadd) = typeadd(i,k)(3:3)
	    pair(3,iadd) = bptadd(i,k)
	    if(k.lt.3) then
	    prtype(iadd)='T'
	    else
	    prtype(iadd)='B'
c	    write(6,*) 'Bifurcated base pair betn.',i,nbpadd(i,k)
	    endif
	    npadd1(npradd) = ibapr(i*2-1)
	    npadd2(npradd) = nbpadd(i,k)
	    if(npadd1(npradd).ne.0) then
	    badd1(npradd) = secnm(ibapr(i*2-1))
	    badd2(npradd) = secnm(nbpadd(i,k))
c	    write(6,72) ibapr(i*2-1),badd1(npradd),nbpadd(i,1),badd2(npradd),
c     1 npradd
	    call getblprm(ibapr(i*2-1),nbpadd(i,k),typeadd(i,k),
     1 bptadd(i,k), bptilt(iadd),bproll(iadd),bptwst(iadd),
     2 bpshft(iadd),bpslid(iadd),bprise(iadd),stk,ybs)
c	    write(8,72) ibapr(i*2-1),badd1(npradd),nbpadd(i,k),badd2(npradd),
c     1 npradd
	    endif
	  endif
 601	continue
 72	format(i5,1x,a1,i5,1x,a1,i5)

C
      IF(NTIM .GT. 1)GO TO 500	
        IF(NGLOB .EQ. 1) THEN
c         WRITE(6,99)
         WRITE(8,99)
 99   FORMAT(/,'Should the molecule be reoriented about a LINEAR GLOBAL
     1Helix Axis? [Y]',$)
c        READ(5,'(A1)')ANSrnt
        IF(ANSrnt.NE.'n'.AND.ANSrnt.NE.'N') THEN
	  ANSrnt = 'Y'
	ELSE
	  ANSrnt = 'N'
	END IF
        WRITE(8,'(5X,A1)')ANSrnt
        IF(ANSrnt.NE.'Y') GO TO 500
C
        END IF
C
           NPTS = NSTEPS
         IF(NGLOB .EQ. 1)THEN
        REWIND 9
	      write(8,89)
              WRITE(8,'(5X,A4)') anorin
        CALL REORNT(IBAPR,CENTRS,NSTEPS,ANORIN,NPOINT,ANSDB)
	 IF(ANORIN.EQ.'O   '.OR.ANORIN.EQ.'C   ') THEN
           CALL LINEFIT(DC,NPOINT,CENTRS)
C
	 ELSE IF(ANSDB.NE.'Y') THEN
           CALL LINASN(DC,NPOINT,CENTRS)
	 ELSE 
	    CALL LINATM(DC,NPOINT,CENTRS)
	 END IF
           ZGX = 0.0D0
           ZGY = 0.0D0
           ZGZ = 1.0D0
           CALL CROSS(DC(1),DC(2),DC(3),ZGX,ZGY,ZGZ,DCR1,DCR2,DCR3)
           ANGL = DACOS(DC(3))
           ANGLE = ANGL * CONV
           DCRS1 = real(DCR1)
           DCRS2 = real(DCR2)
           DCRS3 = real(DCR3)
           CALL ROTMAT(DCRS1,DCRS2,DCRS3,ANGLE,RMAT)
            DO 2006 KJ=1,NAT
              CALL MATMUL(RMAT,CR(KJ,1),CR(KJ,2),CR(KJ,3),TMP)
                DO 2007 KL=1,3
                     CR(KJ,KL) = TMP(KL)
 2007	         CONTINUE
 2006	     CONTINUE
C
           ELSE IF(NGLOB .EQ. 2)THEN
             REWIND 9
             CALL REORNT(IBAPR,CENTRS,NPTS,ANORIN,NPOINT,ANSDB)
             IF(ANORIN.EQ.'O   '.OR.ANORIN.EQ.'C   ') THEN
                CALL  POINT(NPOINT,DX1,DY1,CENTRS)
             ELSE IF(ANSDB.EQ.'Y') THEN
                 CALL PONTDB(NPOINT,CENTRS,DX1,DY1)
             ELSE 
                 CALL PONTSN(NPOINT,CENTRS,DX1,DY1)
             END IF
             DO 2008 KJ=1,NAT
                CR(KJ,1) = CR(KJ,1) - DX1
                CR(KJ,2) = CR(KJ,2) - DY1
 2008        CONTINUE
          END IF
C
 2005   CONTINUE
 500  CONTINUE

      if((ansovl.ne.'N').and.(ansovl.ne.'n'))then
        call overlap(nb,ibapr,ovrlp)
      endif

      CALL WRITING(filecr,TITLE1,INPUT,NTIM,NRUN,ANORIN,filprt,hlxfl)
C-----------------------------------------------------------------------
      IF(NTIM .GT. 1)GO TO 701
      IF(ANSrnt .EQ. 'Y' ) CALL REWRIT(filecr,TITLE1,KKTI)
C   Reoriented atomic coordinates written in xxx.COOR
C-----------------------------------------------------------------------
        IF(NTIM.EQ.1.AND.ANSrnt.EQ.'Y') WRITE(9,1101)
 701    CONTINUE
        IF(NTIM.EQ.2.AND.ANSrnt.EQ.'Y') WRITE(9,1102)
 1101  FORMAT(//,30X,'First strand axes etc.'/20X,3('--------------'))
 1102  FORMAT(//,30X,'Second strand axes etc.'/20X,3('-------------'))
 700  CONTINUE
C-----------------------------------------------------------------------
c      WRITE(6,9993)
      WRITE(8,9993)
9993  FORMAT(/' Are P--P distances and Cylindrical Polar Coordinates req
     1uired? ' /'Only meaningful for reoriented molecule. [Y]',$)
c      READ(5,'(A1)')ANSPP
      IF(ANSPP.NE.'n'.AND.ANSPP.NE.'N') ANSPP = 'Y'
      WRITE(8,'(5X,A1)') ANSPP
      IF(ANSPP.EQ.'Y') THEN
	 ATTYPE = 'P   '
 	 CALL RTHEPH(ANSDB,ATTYPE,nb,ibapr)
	 ATTYPE = 'C1'' '
	 CALL RTHEPH(ANSDB,ATTYPE,nb,ibapr)
      END IF
C	WRITE(6,9996)
C	WRITE(8,9996)
C 9996  FORMAT(/'Are C1''-C1'' distances and Cylindrical Polar Coordinat
C     1es required? [Y]',$)
C	READ(5,'(A1)') ANSCCYL
C	WRITE(8,'(5X,A1)') ANSCCYL
C	IF(ANSCCYL.NE.'N'.AND.ANSCCYL.NE.'n') THEN
C	   ATTYPE = 'C1'' '
C	END IF
C-----------------------------------------------------------------------
c       WRITE(6,9995)
       WRITE(8,9995)
9995  FORMAT(/'Is the torsion angle calculation required? [Y]',$)
c       READ(5,'(A1)') ANSTOR
       IF(ANSTOR.NE.'n'.AND.ANSTOR.NE.'N') ANSTOR = 'Y'
       WRITE(8,'(5X,A1)') ANSTOR
       IF(ANSTOR.EQ.'Y') CALL TORSION(ANSDB,ANSBPN,filprt)
        close(unit=11)
c	close(unit=8)
c	close(unit=9)
	close(unit=2)
       return
C-----------------------------------------------------------------------
	stop
131	write(6,*) 'Unable to open Coordinate file: ',filecr
	stop
132	write(6,*) 'Unable to open Base Pair Information file: ',INPUT
        stop
       END
C-----------------------------------------------------------------------
      SUBROUTINE GETCOR(IUNIT,IERR,ANSPDB)
	parameter ( natm = 999999 )
	parameter ( nrs = 9999 )
	include 'coordinates.h'
	common /nameorg/secnmo(nrs)
	common /nameno/navar,ngvar,ncvar,nuvar
	common /namevar/adevar(200),guavar(200),cytvar(200),uravar(200)
C..
C     THIS ROUTINE IS TAKEN FROM THE PROGRAM CHEM
C     AND MODIFIED
C..   READING COORDINATES FROM IUNIT IN PDB FORMAT
C
C
	external fraction
      CHARACTER*1 ANSPDB
	character*3 resn,adevar,cytvar,guavar,uravar,secnmo
      CHARACTER*80 LINE
      CHARACTER*40 FORMA
      INTEGER RESNO
C
      DATA MAXSEG,MAXAT/nrs,natm/
 5    FORMAT(A80)
 15   FORMAT(12X,A4,1X,A3,2X,I4,4X,3F8.3)
17	format(a3)
C
  6   FORMAT(' IO ERROR READING *PDB* FILE')
 16   FORMAT(' SYSTEM TOO LARGE')
 26   FORMAT(' BAD ORDERING OF RESIDUES')
 36   FORMAT(' NO ATOMS DEFINED WITHIN *PDB* FILE')
C
        open(unit=14,file='/usr/local/bin/AdeVariants.name')
	i=1
	do while(i.le.200) 
	   read(14,17,end=115)adevar(i)
	   i=i+1
	enddo
115	navar=i-1
	i=1
	close(unit=14)
        open(unit=14,file='/usr/local/bin/GuaVariants.name')
	do while(i.le.200)
	   read(14,17,end=116) guavar(i)
	   i=i+1
	enddo
116	ngvar=i-1
	i=1
	close(unit=14)
        open(unit=14,file='/usr/local/bin/CytVariants.name')
	do while(i.le.200)
	   read(14,17,end=117) cytvar(i)
	   i=i+1
	enddo
117	ncvar=i-1
	i=1
	close(unit=14)
        open(unit=14,file='/usr/local/bin/UraVariants.name')
	do while(i.le.200)
	   read(14,17,end=118) uravar(i)
	   i=i+1
	enddo
118	nuvar=i-1
	close(unit=14)

      NW = 6
      REWIND IUNIT
C
      IERR = 0
      NAT = 0
      NSEG = 0
      LRESNO = -1
c        WRITE(6,996)
        WRITE(8,996)
 996  FORMAT(/'Is the data file in the "BROOKHAVEN PDB" format? [Y]',$)
c        READ(5,'(A1)') ANSPDB
	 IF(ANSPDB.NE.'n'.AND.ANSPDB.NE.'N') THEN
	   ANSPDB = 'Y'
	 ELSE
	   ANSPDB = 'N'
	 END IF
        WRITE(8,'(5X,A1)') ANSPDB
        IF(ANSPDB.EQ.'Y') THEN
C
C     LOOP TO READ THE COORDINATES IN PDB FORMAT
C
 10   CONTINUE
          READ(4,5,END=61,ERR=21) LINE
          IF (LINE(1:4).EQ.'ATOM') THEN
            NAT = NAT + 1
            IF (NAT.GT.MAXAT) GO TO 31
              READ(LINE,15) ATNM(NAT),RESN,RESNO,(CR(NAT,L),L=1,3)
              do while(ATNM(NAT)(1:1).EQ.' ') 
                ATNM(NAT)(1:1)=ATNM(NAT)(2:2)
                ATNM(NAT)(2:2)=ATNM(NAT)(3:3)
                ATNM(NAT)(3:3)=ATNM(NAT)(4:4)
                ATNM(NAT)(4:4)=' '
              enddo
	      do while(resn(1:1).eq.' ')
	        resn(1:1)=resn(2:2)
	        resn(2:2)=resn(3:3)
	        resn(3:3)=' '
	      enddo
	      namfound=0
	      nt=1
	      do k=1,navar
	        if(resn.eq.adevar(k)) namfound=1
	      enddo
	      do k=1,ngvar
	        if(resn.eq.guavar(k)) namfound=1
	      enddo
	      do k=1,ncvar
	        if(resn.eq.cytvar(k)) namfound=1
	      enddo
	      do k=1,nuvar
	        if(resn.eq.uravar(k)) namfound=1
                if(resn.eq.'THY'.or.resn.eq.'DT '.or.resn.eq.'T  ')nt=1
	      enddo
              if(resn.eq.'QUO') namfound=1
              if(resn.eq.'PSU') namfound=1
              if(resn.eq.'FHU') namfound=1
              if(resn.eq.'N6G') namfound=1

	      if(namfound.eq.1) then
	       
	        IRES(NAT) = RESNO
                IF (RESNO.NE.LRESNO) THEN
                  NSEG = NSEG + 1
                  IF (NSEG.GT.MAXSEG) GO TO 31
                    SECSQ(NSEG) = NSEG
                    SECNM(NSEG) = RESN
		    if(resn.eq.'THY'.or.resn.eq.'DT '.or.resn.eq.'T  ')
     1    secnmo(nseg)='THY'
                    ISTS(NSEG) = NAT
                    LRESNO = RESNO
                ENDIF
                IENS(NSEG) = NAT
	      else
	       nat=nat-1
	      endif
          ENDIF
      GO TO 10
C     END OF THE READING LOOP
 61   CONTINUE
        ELSE
            WRITE(6,997)
            WRITE(8,997)
 997  FORMAT(/' If not, type in the input data format'/' Field specific
     1ations required: Atom Name, Res.Name, Res.No., XYZ Coords.')
            READ(5,'(A40)') FORMA
            WRITE(8,'(A40)') FORMA
            WRITE(6,998)
            WRITE(8,998)
 998  FORMAT(/' TYPE UNIT CELL PARAMETERS: A,B,C AND ALPHA,BETA,GAMA')
C
            READ(5,*) AC,BC,CC,ALP,BET,GAM
            WRITE(8,911) AC,BC,CC,ALP,BET,GAM
 911  FORMAT(/2X,'UNIT CELL PARAMETERS :: A =',F8.2,' B =',F8.2,' C =',
     1       F8.2,/21X,' APLHA =',F8.2,' BETA =',F8.2,' GAMA =',F8.2)
 101    CONTINUE
              READ(4,5,END=161,ERR=21) LINE
C
              NAT = NAT + 1
              IF (NAT.GT.MAXAT) GO TO 31
              READ(LINE,FORMA) ATNM(NAT),RESN,RESNO,(CR(NAT,L),L=1,3)
              IF(ATNM(NAT)(1:1).EQ.' ') THEN
                 ATNM(NAT)(1:1)=ATNM(NAT)(2:2)
                 ATNM(NAT)(2:2)=ATNM(NAT)(3:3)
                 ATNM(NAT)(3:3)=ATNM(NAT)(4:4)
                 ATNM(NAT)(4:4)=' '
              ENDIF
	      IRES(NAT) = RESNO
C
              IF (RESNO.NE.LRESNO) THEN
                  NSEG = NSEG + 1
                  IF (NSEG.GT.MAXSEG) GO TO 31
                  SECSQ(NSEG) = NSEG
                  SECNM(NSEG) = RESN
                  ISTS(NSEG) = NAT
                  LRESNO = RESNO
              ENDIF
C
              IENS(NSEG) = NAT
      GO TO 101
C     END OF THE READING LOOP
 161  CONTINUE
        CALL FRACTION(AC,BC,CC,ALP,BET,GAM)
        END IF
      IF (NAT.EQ.0) GO TO 51
      IERR = 0
      RETURN
C
 21   CONTINUE
      WRITE (NW,6)
      IERR = 1
      RETURN
C
 31   CONTINUE
      WRITE (NW,16)
      IERR = 1
      RETURN
C
C 41   CONTINUE
C     WRITE (NW,26)
C
 51   CONTINUE
      WRITE (NW,36)
      IERR = 1
      RETURN
      END
 
      SUBROUTINE PICATP(IBAPR1,IBAPR2,IPAIR,IBASE1,IBASE2,NB1,
     1                  NB2,INOHP,BASE1,BASE2)
	parameter ( natm = 999999 )
	parameter ( nrs = 9999 )
C
	include 'coordinates.h'
C
      DIMENSION IPAIR(4),ATOMPA(9)
	common /nameorg/secnmo(nrs)
	common /nameno/navar,ngvar,ncvar,nuvar
	common /namevar/adevar(200),guavar(200),cytvar(200),uravar(200)
      DIMENSION IBASE1(10),IBASE2(10)
C
      CHARACTER*1 PRIME
      CHARACTER*3 TMP,adevar,cytvar,guavar,uravar,secnmo
      CHARACTER*3 BASE1,BASE2
      CHARACTER*4 ATOMPA
C
      DATA ATOMPA/'C4'' ','C3'' ','O4'' ','C2'' ','C1'' ','N9  ','N1  ',
     1'C8  ','C6  '/
      DATA PRIME/''''/
C
      DO 2000 KL=1,4
      IPAIR(KL) = 0
 2000	CONTINUE
C
C     TO CHANGE THE C1* etc to C1' AND O1' TO O4' TO MATCH
C     THE ATOM NAMES FROM PDB AND AMBER
C
      DO 345 IA = 1, NAT
C      IF(ATNM(IA)(1:3).EQ.'O1'' ')ATNM(IA)=ATOMPA(3)
      IF(ATNM(IA)(3:3).EQ.'*')ATNM(IA)(3:3)=PRIME
 345  CONTINUE
 555  FORMAT(3X,I5,30(2X,A1))
c	write(*,*) 'No. of Gua, Ade, Cyt and Ura variants',ngvar,navar,
c     1 ncvar,nuvar
      DO 445 IR = 1,NSEG
	do i=1,ngvar
	   if(secnm(ir).eq.guavar(i)) secnm(ir)='GUA'
	enddo
	do i=1,navar
	   if(secnm(ir).eq.adevar(i)) secnm(ir)='ADE'
	enddo
	do i=1,ncvar
	   if(secnm(ir).eq.cytvar(i)) secnm(ir)='CYT'
	enddo
	do i=7,nuvar
	   if(secnm(ir).eq.uravar(i)) secnm(ir)='URA'
	enddo
	do i=1,6
         if(secnm(ir).eq.uravar(i).and.secnmo(ir).eq.'THY')
     1    secnm(ir)='THY'
	enddo
	if(secnm(ir).eq.'PSU') secnm(ir)='URA'
	if(secnm(ir).eq.'FHU') secnm(ir)='URA'
	if(secnm(ir).eq.'N6G') secnm(ir)='GUA'
	if(secnm(ir).eq.'QUO') secnm(ir)='GUA'
	if((secnm(ir).ne.'ADE').and.(secnm(ir).ne.'GUA').and.
     1  (secnm(ir).ne.'CYT').and.(secnm(ir).ne.'URA').and.
     1  (secnm(ir).ne.'THY'))then
	write(*,337)ir,secnm(ir)
	endif
337	format('MODRES',3X,I5,5x,A3,4x,A3)

c	write(*,*) ir,secnm(ir) 
445   CONTINUE
      ITYPE1 = 2
      ITYPE2 = 2
      IF(SECNM(IBAPR1).EQ.'ADE'.OR.SECNM(IBAPR1).EQ.'GUA') ITYPE1=1
      IF(SECNM(IBAPR2).EQ.'ADE'.OR.SECNM(IBAPR2).EQ.'GUA') ITYPE2=1
      ISTR = ISTS(IBAPR1)
      BASE1 = SECNM(IBAPR1)
      IEND = IENS(IBAPR1)
      DO 100 IAT = ISTR,IEND
        IF(ITYPE1.EQ.1) THEN
          IF(ATNM(IAT)(1:2).EQ.ATOMPA(6)(1:2)) IPAIR(2)=IAT
          IF(ATNM(IAT)(1:2).EQ.ATOMPA(8)(1:2)) IPAIR(1)=IAT
        ELSE IF(ITYPE1.EQ.2) THEN
          IF(ATNM(IAT)(1:2).EQ.ATOMPA(7)(1:2)) IPAIR(2)=IAT
          IF(ATNM(IAT)(1:2).EQ.ATOMPA(9)(1:2)) IPAIR(1)=IAT
        ENDIF
100   CONTINUE
C
c      write(*,*)istr,iend,itype1
      CALL BASPIC(ISTR,IEND,ITYPE1,IBASE1,NB1)
C
      IF(IBAPR2.NE.0) THEN
C
      ISTR = ISTS(IBAPR2)
      IEND = IENS(IBAPR2)
      BASE2 = SECNM(IBAPR2)
      DO 101 IAT = ISTR,IEND
        IF(ITYPE2.EQ.1) THEN
          IF(ATNM(IAT)(1:2).EQ.ATOMPA(6)(1:2)) IPAIR(4)=IAT
          IF(ATNM(IAT)(1:2).EQ.ATOMPA(8)(1:2)) IPAIR(3)=IAT
        ELSE IF(ITYPE2.EQ.2) THEN
          IF(ATNM(IAT)(1:2).EQ.ATOMPA(7)(1:2)) IPAIR(4)=IAT
          IF(ATNM(IAT)(1:2).EQ.ATOMPA(9)(1:2)) IPAIR(3)=IAT
        ENDIF
101   CONTINUE
C
      CALL BASPIC(ISTR,IEND,ITYPE2,IBASE2,NB2)
C
      ELSE
C
      IPAIR(3) = 0
      BASE2 = '   '
C
       DO 2001 IAT=ISTR,IEND
          IF(ATNM(IAT)(1:3).EQ.ATOMPA(5)(1:3)) IPAIR(4) = IAT
 2001	CONTINUE
      ENDIF
C
      RETURN
      END
 
C----------------------------------------------------------------------
      SUBROUTINE BASPIC(ISTR,IEND,ITYPE,IATPIC,NB)
	parameter ( natm = 999999 )
	parameter ( nrs = 9999 )
	include 'coordinates.h'
C
      DIMENSION ATPUR(10),ATPYR(7),IATPIC(10)
      CHARACTER*4 ATPUR,ATPYR
 
      DATA ATPUR/'C1'' ','N9  ','C8  ','N7  ','C5  ','C6  ','N1  ',
     1            'C2  ','N3  ','C4  '/
      DATA ATPYR/'C1'' ','N1  ','C6  ','C5  ','C4  ','N3  ','C2  '/
      DATA NPUR,NPYR/10,7/
C
C     WRITE(9,887) ITYPE
887   FORMAT(' ITYPE ' , I4/)
      GO TO (20,30),ITYPE
C
C     PURINE  BASE ATOMS
C
20    CONTINUE
      NB=NPUR
      DO 100 IA = ISTR,IEND
      DO 100 IP = 1,NPUR
      IF(ATNM(IA)(1:3).EQ.ATPUR(IP)(1:3)) IATPIC(IP)=IA
100   CONTINUE
     
C     WRITE(9,57)(IATPIC(IC),IC=1,NPUR)
57    FORMAT(20I4)

      RETURN
30    CONTINUE
      NB=NPYR
      DO 200 IA = ISTR,IEND
      DO 200 IP = 1,NPYR
      IF(ATNM(IA)(1:3).EQ.ATPYR(IP)(1:3)) IATPIC(IP)=IA
200   CONTINUE
C     WRITE(9,57)(IATPIC(IC),IC=1,NPYR)
      RETURN
      END
 
C----------------------------------------------------------------------!
      SUBROUTINE FINDVEC(IPS,INS,NAT1,NAT2,NAT3,NAT4,NOBASE,ANSOR,ANSC1,
     1ANSrnt,NTIM,ibapr1,ibapr3,ybs1,ybs2,type1,type2,type3,type4,
     2cystrns,cystrns2)
C----------------------------------------------------------------------!
	parameter ( natm = 999999 )
	parameter ( nrs = 9999 )
        include 'coordinates.h'
        include 'parameters.h'
      DOUBLE PRECISION EL,EM,EN,EL2,EM2,EN2,XAX,XAXIS,YAXIS,DOTPR,BM,
     1   BMN,WGLT,WGLL,WGTILT,WGROLL,PROJ,ELT,EMT,ENT,ELN,EMN,ENN,ZAXIS,
     2   BPDIR,BPNORML,B,C8C6,YC8C6, EL1,EM1,EN1,SM,WTILT,WROLL,WGLT1,
     3   WGRL1,BNORM,EK1,EK2,EK3,VEC,TEMPX,TEMPY,xbase1,xbase3,ybase1,
     4   ybase3,zbase1,zbase3,BPNORMLX,BPNORMLY,MBPNORML,AMFO,ybs1,ybs2,
     5   possblz,px,py,pz
      DIMENSION X(3,15),ATOM(15),NSP(15),ORIG(3),PROJ(3),XYZ1(3),
     1   XYZ2(3),XYZ3(3),XYZ4(3),DVA(15),N(15),PROJ1(2),B(2,3),AX(3),
     2   IPS(4),C(3),XAX(3),XAXIS(2,3),YAXIS(2,3),BMN(3),INS(4),
     3   NOATM(2),VEC(3),X1(3,15),ATOM1(15),NSP1(15),N1(15),XC8(3),
     4   XC6(3),RC8(3),RC6(3),SM(3),BM(3),BNORM(3,4),ORIGIN(2,3),
     5   TEMPX(3),TEMPY(3),R(3,3),XC1P(3),RC1P(3),C8C1(3),C6C1(3),
     6   XC1P2(3),ZAXIS(2,3),BPDIR(3,4),ORIG1(3),ORIG2(3),XC81(3),
     7   XC82(3),XC1P1(3),C8C6(3),YC8C6(3),XMID(3),orig3(3),xbase1(3),
     8   xbase3(3),ybase1(3),ybase3(3),zbase1(3),zbase3(3),ybs1(2,3),
     9   ybs2(2,3),possblZ(3)
      COMMON /BPZAXS/BPNORML(nrs,3),BPNORMLX(nrs,3),BPNORMLY(nrs,3),
     1     MBPNORML(nrs,3),AMFO(nrs,3)
      COMMON /BASINF/BASE(nrs),pair(3,nrs),prtype(nrs)
      COMMON /ATOMS/IBASE1(10),IBASE2(10),IBASE3(10),IBASE4(10),NB,
     1 ibapr(nrs)
      LOGICAL DOUBL1,DOUBL2,STKINV,RETCAL
      CHARACTER*4 ATOM,ATOM1
      CHARACTER*3 BASE
      CHARACTER*1 TMPBAS,ANSOR,ANSC1,ANSrnt,pair,prtype,type1,type2,
     1type3,type4,cystrns,cystrns2
C
      STKINV = .FALSE.
C
	write(9,*) '==================================================='
	write(9,*)'Calculating parameters betn. IBAPR1 & IBAPR3',ibapr1
     1 ,ibapr3
      NRTCAL = 0
      CONV=180.0/3.141592654
      DOUBL1 = .TRUE.
      DOUBL2 = .TRUE.
C
      DO 2000 J=1,NAT1
        IF(IBASE1(J).EQ.INS(1)) NC8 = J
        IF(IBASE1(J).EQ.INS(2)) NC6 = J
        IF(IBASE1(J).EQ.IPS(1)) NN9 = J
        IF(IBASE1(J).EQ.IPS(2)) NN1 = J
        KJ = IBASE1(J)
        ATOM(J) = ATNM(KJ)
        N(J) = 0
        DO 2001 KC =1,3
           X(KC,J) = CR(KJ,KC)
 2001   CONTINUE
C
        IF(ATOM(J).EQ.'C1'' ') THEN
           DO 2002 KC=1,3
              XC1P(KC) = X(KC,J)
              XC1P1(KC) = X(KC,J)
 2002      CONTINUE
         ELSE
           N(J) = J
         ENDIF
C
 2000  CONTINUE
C
       DO 2003 KC = 1,3
           XC8(KC) = CR(INS(1),KC)
           XC81(KC) = CR(INS(1),KC)
           XC6(KC) = CR(INS(2),KC)
           XYZ1(KC) = CR(IPS(1),KC)
           XYZ4(KC) = CR(IPS(2),KC)
           XYZ2(KC) = CR(INS(1),KC)
           XYZ3(KC) = CR(INS(2),KC)
 2003  CONTINUE
C 
      IF(INS(2).NE.0) THEN
        DO 2004 I=1,NAT2
          IF(IBASE2(I).EQ.INS(1)) NC8 = I
          IF(IBASE2(I).EQ.INS(2)) NC6 = I
          IF(IBASE2(I).EQ.IPS(1)) NN9 = I
          IF(IBASE2(I).EQ.IPS(2)) NN1 = I
          KJ = IBASE2(I)
          ATOM1(I) = ATNM(KJ)
          N1(I) = 0
          DO 2005 KC=1,3
            X1(KC,I) = CR(KJ,KC)
2005      CONTINUE
          IF(ATOM1(I).EQ.'C1'' ') THEN
            DO 2006 KC=1,3
              RC1P(KC) = X1(KC,I)
2006        CONTINUE
          ELSE
            N1(I) = I
          ENDIF
2004	CONTINUE
C --------------------------------------------------------------------
C   ASSIGNED THE ATOM NOS. FOR FIXING Y-AXIS .
C   NOW FINDING THE COORDINATES OF BASE PAIR ORIGIN AND DIRECTION
C   RATIOS OF THE "Y"-AXIS.
C --------------------------------------------------------------------
        SUMANG = 0.0
        SUMDS1 = 0.0
        SUMDS2 = 0.0
        DO 2008 KJ=1,3
          C8C1(KJ) = XC1P(KJ) - XC8(KJ)
          C6C1(KJ) = RC1P(KJ) - XC6(KJ)
          SUMDS1 = C8C1(KJ) * C8C1(KJ) + SUMDS1
          SUMDS2 = C6C1(KJ) * C6C1(KJ) + SUMDS2
2008	CONTINUE
        SUMDS1 = SQRT(SUMDS1)
        SUMDS2 = SQRT(SUMDS2)
        DO 2009 KJ=1,3
          C8C1(KJ) = C8C1(KJ)/SUMDS1
          C6C1(KJ) = C6C1(KJ)/SUMDS2
          BPDIR(KJ,1) = C8C1(KJ) * 1.0
          BPDIR(KJ,2) = C6C1(KJ) * 1.0
2009	CONTINUE
	XN9C1 = 0.0
	XN1C1 = 0.0
	C1C1DS = 0.0
        C8C6DIS = 0.0
	DSN1C1 = 0.0
	DSN9C1 = 0.0
        DO 2010 KJ=1,3
          C8C6DIS = C8C6DIS + (XC8(KJ) - XC6(KJ)) * (XC8(KJ) -XC6(KJ))
	  C1C1DS = C1C1DS + (XC1P(KJ) - RC1P(KJ)) * (XC1P(KJ) -RC1P(KJ))
	  XN9C1 = (X(KJ,NN9) - XC1P(KJ)) * (RC1P(KJ)-XC1P(KJ))+ XN9C1
	  XN1C1 = (X1(KJ,NN1) - RC1P(KJ)) * (XC1P(KJ)-RC1P(KJ))+ XN1C1
	DSN9C1 = DSN9C1 + (X(KJ,NN9) - XC1P(KJ)) * (X(KJ,NN9) -XC1P(KJ))
	DSN1C1 = DSN1C1 + (X1(KJ,NN1) - RC1P(KJ)) *(X1(KJ,NN1)-RC1P(KJ))
2010    CONTINUE
	DSN9C1 = SQRT(DSN9C1)
	DSN1C1 = SQRT(DSN1C1)
	OPENC1(LNDEX + 1) = SQRT(C1C1DS)
        OPENDS(LNDEX + 1) = SQRT(C8C6DIS)
	ANGBN9(LNDEX+1) = ACOS(XN9C1/(DSN9C1 * OPENC1(LNDEX+1))) * CONV
	ANGBN1(LNDEX+1) = ACOS(XN1C1/(DSN1C1 * OPENC1(LNDEX+1))) * CONV
        DO 2011 KJ=1,3
          IF(ANSOR.EQ.'Y') THEN
            ORIG1(KJ)=0.0
            ORIG2(KJ)=0.0
            DO 2012 NN=1,NAT1 
              IF(ATOM(NN) .NE. 'C1'' ')ORIG1(KJ) = ORIG1(KJ) + X(KJ,NN)
2012	    CONTINUE
            ORIG1(KJ) = ORIG1(KJ)/(NAT1-1)
            DO 2013 NN=1,NAT2
              IF(ATOM1(NN).NE.'C1'' ')ORIG2(KJ) = ORIG2(KJ) + X1(KJ,NN)
2013	    CONTINUE
            ORIG2(KJ) = ORIG2(KJ)/(NAT2-1)
            ORIG (KJ) = (ORIG1(KJ) + ORIG2(KJ)) * 0.5
          ELSE
            ORIG(KJ) = (XC8(KJ) + XC6(KJ)) * 0.5
          END IF
          PROJ(KJ) = XC8(KJ) - XC6(KJ)
2011	CONTINUE
        IF(ANSC1.EQ.'Y') THEN
          DO 2014 KKK=1,3
            PROJ(KKK) = XC1P(KKK) - RC1P(KKK)
	    XMID(KKK) = (XC1P(KKK) + RC1P(KKK)) * 0.5
	    C8C6(KKK) = XC8(KKK) - XC6(KKK)
 2014	  CONTINUE
	  CALL NORMAL(C8C6(1),C8C6(2),C8C6(3),YC8C6(1),YC8C6(2),
     1    YC8C6(3))
          CALL NORMAL(PROJ(1),PROJ(2),PROJ(3),YAXIS(1,1),YAXIS(1,2),
     1    YAXIS(1,3))
	  CONSTN=(YAXIS(1,1)*(XMID(1)-XC8(1))+YAXIS(1,2)*(XMID(2)-XC8(2))
     1    +YAXIS(1,3)*(XMID(3)-XC8(3)))/(YAXIS(1,1)*YC8C6(1)+YAXIS(1,2)*
     2    YC8C6(2) + YAXIS(1,3)*YC8C6(3))
	  XYZMD1 = XC8(1) + YC8C6(1) * CONSTN
	  XYZMD2 = XC8(2) + YC8C6(2) * CONSTN
	  XYZMD3 = XC8(3) + YC8C6(3) * CONSTN
	ELSE
          CALL NORMAL(PROJ(1),PROJ(2),PROJ(3),YAXIS(1,1),YAXIS(1,2),
     1    YAXIS(1,3))
        END IF
	IF((ANSOR.EQ.'Y').OR.(ANSC1.NE.'Y')) THEN
 	  ORIGIN(1,1) = ORIG(1)
          ORIGIN(1,2) = ORIG(2)
	  ORIGIN(1,3) = ORIG(3)
	ELSE
	  ORIGIN(1,1) = XYZMD1
	  ORIGIN(1,2) = XYZMD2
	  ORIGIN(1,3) = XYZMD3
	END IF
        WRITE(9,'(''   BASE-PAIR ORIGIN-1'',3F10.4)')
     1  (ORIGIN(1,KJ),KJ=1,3)
        WRITE(21,'(''BASE-PAIR ORIGIN-1'',3F10.4)') 
     1  (ORIGIN(1,KJ),KJ=1,3)

C----------------------------------------------------------------------
C FIND L,M,N, OF THE TWO BASE NORMALS & B.P. NORMAL IN FIRST base-pair
C----------------------------------------------------------------------

	if(nat1.ge.5) then
          CALL PLANEB(X,NAT1,N,1,DVA,EL,EM,EN,DETR,SDV1)
          BNORM(1,1) = EL
          BNORM(2,1) = EM
          BNORM(3,1) = EN
	else
	  write(6,801) nat1
 801	format('BAD Selection of BASE, only',I3,' atoms found')
	  stop
	endif
C
	if(nat2.ge.5) then
          CALL PLANEB(X1,NAT2,N1,1,DVA,EL1,EM1,EN1,DETR,SDV2)
	else
	  write(6,801) nat2
	  stop
	endif
        ANG = EL * EL1 + EM * EM1 + EN * EN1
        IF(ANG.LT.0.0) THEN
          EL1 = -EL1
          EM1 = -EM1
          EN1 = -EN1
        ENDIF
        BNORM(1,2) = EL1
        BNORM(2,2) = EM1
        BNORM(3,2) = EN1

C MEAN OF THE BASE-NORMALS TO GET B.P. NORMAL
        EL1 = EL + EL1
        EM1 = EM + EM1
        EN1 = EN + EN1
        CALL NORMAL(EL1,EM1,EN1,EL,EM,EN)
      ELSE
        IF(NC8.EQ.0.AND.NC6.NE.0) THEN
           NC8 = NC6
           DO 2015 KL = 1,3
             XC8(KL) = XC6(KL)
2015       CONTINUE
        ENDIF
c       CALL SINGAXS(IPS,N,NAT1,NC8,X,1,XAXIS,YAXIS,EL,EM,EN,IBASE1,SD1)
	call findybs(ibapr1,xbase1,ybase1,zbase1,ang1,orig1,'W','C')
	do i=1,3
	  xaxis(1,i) = ybase1(i)
	  yaxis(1,i) = zbase1(i)
	  origin(1,i) = orig1(i)
	enddo
	el = xbase1(1)
	em = xbase1(2)
	en = xbase1(3)
        DOUBL1 = .FALSE.
      ENDIF
C----------------------------------------------------------------------
C  SECOND BASE/base-pair OF THE DOUBLET
C----------------------------------------------------------------------
      DO 2016 J=1,NAT3
        IF(IBASE3(J).EQ.INS(3)) NC8 = J
        IF(IBASE3(J).EQ.INS(4)) NC6 = J
        IF(IBASE3(J).EQ.IPS(3)) NN9 = J
        IF(IBASE3(J).EQ.IPS(4)) NN1 = J
        KJ = IBASE3(J)
        ATOM(J) = ATNM(KJ)
        N(J) = 0
        DO 2017 KC = 1,3
          X(KC,J) = CR(KJ,KC)
2017    CONTINUE
        IF(ATOM(J).EQ.'C1'' ') THEN
          DO 2018 KC=1,3
            XC1P(KC) = X(KC,J)
            XC1P2(KC) = X(KC,J)
            possblz(kc)=xc1p2(kc)-xc1P1(kc)
2018      CONTINUE
          call normal(possblz(1),possblz(2),possblz(3),px,py,pz)
          possblz(1)=px
          possblz(2)=py
          possblz(3)=pz
        ELSE
          N(J) = J
        END IF
2016  CONTINUE
      write(9,*) 'Possible Z direction from C1''',possblz
      DO 2019 KC = 1,3
        RC8(KC) = CR(INS(3),KC)
        XC82(KC) = CR(INS(3),KC)
        RC6(KC) = CR(INS(4),KC)
        XYZ1(KC) = CR(IPS(3),KC)
        XYZ4(KC) = CR(IPS(4),KC)
        XYZ2(KC) = CR(INS(3),KC)
        XYZ3(KC) = CR(INS(4),KC)
2019  CONTINUE
      IF(INS(4).NE.0) THEN
        DO 2020 J=1,NAT4
          IF(IBASE4(J).EQ.INS(3)) NC8 = J
          IF(IBASE4(J).EQ.INS(4)) NC6 = J
          IF(IBASE4(J).EQ.IPS(3)) NN9 = J
          IF(IBASE4(J).EQ.IPS(4)) NN1 = J
          KJ = IBASE4(J)
          ATOM1(J) = ATNM(KJ)
          DO 2021 KC = 1,3
            X1(KC,J) = CR(KJ,KC)
2021      CONTINUE
C         N1(J) = 0
          IF(ATOM1(J).EQ.'C1'' ') THEN
            DO 2022 KC=1,3
              RC1P(KC) = X1(KC,J)
C  SUKANYA - MODIFIED
              N1(J)=0
2022        CONTINUE
          ELSE
            N1(J) = J
          END IF
2020    CONTINUE
        DO 2023 KJ=1,3
          C8C1(KJ) = XC1P(KJ) - RC8(KJ)
          C6C1(KJ) = RC1P(KJ) - RC6(KJ)
2023    CONTINUE
        SUMANG = 0.0
        SUMDS1 = 0.0
        SUMDS2 = 0.0
        DO 2024 KJ=1,3
          SUMDS1 = C8C1(KJ) * C8C1(KJ) + SUMDS1
          SUMDS2 = C6C1(KJ) * C6C1(KJ) + SUMDS2
2024	CONTINUE
        SUMDS1 = SQRT(SUMDS1)
        SUMDS2 = SQRT(SUMDS2)
        DO 2025 KJ=1,3
          C8C1(KJ) = C8C1(KJ)/SUMDS1
          C6C1(KJ) = C6C1(KJ)/SUMDS2
          BPDIR(KJ,3) = C8C1(KJ) * 1.0
          BPDIR(KJ,4) = C6C1(KJ) * 1.0
2025	CONTINUE
	DSN9C1 = 0.0
	DSN1C1 = 0.0
	XN1C1 = 0.0
	XN9C1 = 0.0
	C1C1DS = 0.0
        C8C6DIS = 0.0
        DO 2026 KJ=1,3
          C8C6DIS=C8C6DIS+(RC8(KJ)-RC6(KJ))*(RC8(KJ)-RC6(KJ))
	  XN9C1=XN9C1+(X(KJ,NN9)-XC1P(KJ))*(RC1P(KJ)-XC1P(KJ))
	  DSN9C1=DSN9C1+(X(KJ,NN9)-XC1P(KJ))*(X(KJ,NN9)-XC1P(KJ))
	  XN1C1=XN1C1+(X1(KJ,NN1)-RC1P(KJ))*(XC1P(KJ)-RC1P(KJ))
	  DSN1C1=DSN1C1+(X1(KJ,NN1)-RC1P(KJ))*(X1(KJ,NN1)-RC1P(KJ))
	  C1C1DS=C1C1DS+(XC1P(KJ)-RC1P(KJ))*(XC1P(KJ)-RC1P(KJ))
2026    CONTINUE
	DSN1C1 = SQRT(DSN1C1)
	DSN9C1 = SQRT(DSN9C1)
	OPENC1(LNDEX+2) = SQRT(C1C1DS)
        OPENDS(LNDEX+2) =  SQRT(C8C6DIS)
	ANGBN9(LNDEX+2) = ACOS(XN9C1/(DSN9C1 * OPENC1(LNDEX+2))) * CONV
	ANGBN1(LNDEX+2) = ACOS(XN1C1/(DSN1C1 * OPENC1(LNDEX+2))) * CONV
C----------------------------------------------------------------------
C      CHECKING  WHETHER THE MOLECULE GOES UP OR DOWN IN Z.
C----------------------------------------------------------------------
        DO 2027 K=1,3
          ORIG2(K) = (RC8(K) + RC6(K)) * 0.5
          VEC(K) = ORIG2(K) - ORIG(K)
2027	CONTINUE
      ELSE
        DOUBL2 = .FALSE.
        IF(NC8.EQ.0.AND.NC6.NE.0) THEN
          NC8 = NC6
          DO 2028 KL=1,3
            RC8(KL) = RC6(KL)
2028	  CONTINUE
        ENDIF
        DO 2029 KJ=1,3
          VEC(KJ) = RC8(KJ) - XC8(KJ)
2029    CONTINUE
      ENDIF
      CALL NORMAL(VEC(1),VEC(2),VEC(3),EK1,EK2,EK3)
      ANG = EL * EK1 + EM * EK2 + EN * EK3
      IF(ANG.LT.0.0) THEN
        EL = -EL
        EM = -EM
        EN = -EN
        IF(DOUBL1) THEN
c----------------------------------------------------------------
c   done on 6.july 2005 to account for bases which form cusps, or are
c   not stacked, opppositely oriented, specially seen in RNA
c----------------------------------------------------------------------
          DO 2030 KL=1,3
            BNORM(KL,1) = -BNORM(KL,1)
            BNORM(KL,2) = -BNORM(KL,2)
 2030     CONTINUE
        ELSE
c         DO 2031 KL=1,3
c           YAXIS(1,KL) = -YAXIS(1,KL)
c2031     CONTINUE
        ENDIF
      ENDIF

      IF(.NOT.DOUBL1) THEN
c       THETA = 10.0601 				! bhatta July 2005
c       IF(NAT1.EQ.10) THETA = 12.6921
c       ELS = EL * 1.0
c       EMS = EM * 1.0
c       ENS = EN * 1.0
c       CALL ROTMAT(ELS,EMS,ENS,THETA,R)
c       CALL ROTVEC(R,XAXIS(1,1),XAXIS(1,2),XAXIS(1,3),TEMPX)
c       CALL ROTVEC(R,YAXIS(1,1),YAXIS(1,2),YAXIS(1,3),TEMPY)
c       DO 2032 KJ=1,3				! bhatta July 2005
c         YAXIS(1,KJ) = TEMPY(KJ)
c         XAXIS(1,KJ) = TEMPX(KJ)
c         ORIGIN(1,KJ) = XC8(KJ) - YAXIS(1,KJ) * 4.9122
c2032   CONTINUE
      ELSE
        CALL CROSS(YAXIS(1,1),YAXIS(1,2),YAXIS(1,3),EL,EM,EN,XAX(1),
     1  XAX(2),XAX(3))
        CALL NORMAL(XAX(1),XAX(2),XAX(3),XAXIS(1,1),XAXIS(1,2),
     1  XAXIS(1,3))
         call direcx(xaxis,1,possblz,el,em,en)
c        CALL DIRECX(XAXIS,1,XC1P1,XC81,YAXIS,ybs1,NDX,type1,type2,
c     1  cystrns)    !Sukanya, 12 Nov, 2013
        IF(NDX.EQ.1) THEN
          DO 2033 KJ=1,3
            BNORM(KJ,1) = -BNORM(KJ,1)
            BNORM(KJ,2) = -BNORM(KJ,2)
2033	  CONTINUE
C         TMPBAS = BASE(NOBASE)(1:1)
C         BASE(NOBASE)(1:1) = BASE(NOBASE)(3:3)
C         BASE(NOBASE)(3:3) = TMPBAS
        END IF
      END IF
C-----------------------------------------------------------------------
C   END OF  FIRST BASE/base-pair CALC.
C   START OF DOUBLE-STRANDED SECOND base-pair CALC.
C-----------------------------------------------------------------------
      IF(DOUBL2) THEN
        DO 2034 I=1,3
          ORIG1(I) = 0.0
          ORIG2(I) = 0.0
          IF(ANSOR.EQ.'Y') THEN
            DO 2035 NN=1,NAT3
              IF(ATOM(NN) .NE. 'C1'' ')ORIG1(I) = ORIG1(I) + X(I,NN)
2035	    CONTINUE
            ORIG1(I) = ORIG1(I)/(NAT3-1)
            DO 2036 NN=1,NAT4
               IF(ATOM1(NN) .NE.'C1'' ')ORIG2(I) = ORIG2(I) + X1(I,NN)
2036	    CONTINUE
            ORIG2(I) = ORIG2(I)/(NAT4-1)
            ORIG(I) = (ORIG1(I) + ORIG2(I)) * 0.5
          ELSE
            ORIG(I) = (RC8(I) + RC6(I)) * 0.5
          END IF
          PROJ(I) = RC8(I) - RC6(I)
2034    CONTINUE
        IF(ANSC1.EQ.'Y') THEN
          DO 2037 KKK=1,3
            PROJ(KKK) = XC1P(KKK) - RC1P(KKK)
            XMID(KKK) = (XC1P(KKK) + RC1P(KKK)) * 0.5
	    C8C6(KKK) = RC8(KKK) - RC6(KKK)
 2037	  CONTINUE
	  CALL NORMAL(C8C6(1),C8C6(2),C8C6(3),YC8C6(1),YC8C6(2),
     1    YC8C6(3))
          CALL NORMAL(PROJ(1),PROJ(2),PROJ(3),YAXIS(2,1),YAXIS(2,2),
     1    YAXIS(2,3))
	  CONSTN=(YAXIS(2,1)*(XMID(1)-RC8(1))+YAXIS(2,2)*(XMID(2)
     1    -RC8(2))+YAXIS(2,3)*(XMID(3)-RC8(3)))/(YAXIS(2,1)*YC8C6(1)
     2    + YAXIS(2,2)*YC8C6(2) + YAXIS(2,3)*YC8C6(3))
	  XYZMD1 = RC8(1) + YC8C6(1) * CONSTN
	  XYZMD2 = RC8(2) + YC8C6(2) * CONSTN
	  XYZMD3 = RC8(3) + YC8C6(3) * CONSTN
	ELSE
          CALL NORMAL(PROJ(1),PROJ(2),PROJ(3),YAXIS(2,1),YAXIS(2,2),
     1    YAXIS(2,3))
        END IF
640     FORMAT(3X,'BASE-PAIR ORIGIN-2',3F10.4)
	IF((ANSOR.EQ.'Y').OR.(ANSC1.NE.'Y')) THEN
          ORIGIN(2,1) = ORIG(1)
   	  ORIGIN(2,2) = ORIG(2)
	  ORIGIN(2,3) = ORIG(3)
	ELSE
	  ORIGIN(2,1) = XYZMD1
	  ORIGIN(2,2) = XYZMD2
	  ORIGIN(2,3) = XYZMD3
	END IF
        WRITE(9,640) (ORIGIN(2,KJ),KJ=1,3)
        WRITE(21,'(''BASE-PAIR ORIGIN-2'',3F10.4)') (ORIGIN(2,KJ),
     1 KJ=1,3)
C ----------------------------------------------------------------------
C     FINDING BASE & base-pair NORMALS OF THE TOP PLANE
C ----------------------------------------------------------------------
	if(nat3.ge.5) then
          CALL PLANEB(X,NAT3,N,1,DVA,EL2,EM2,EN2,DETR,SDV3)
	else
	  write(6,801) nat3
 	  stop
	endif
        BNORM(1,3) = EL2
        BNORM(2,3) = EM2
        BNORM(3,3) = EN2
	if(nat4.ge.5) then
C ERRONEOUS ZONE
c	  write(*,*)NAT4
c	  write(*,436)(N1(IPH),IPH=1,NAT4)
c436	  format(I15)
          CALL PLANEB(X1,NAT4,N1,1,DVA,EL1,EM1,EN1,DETR,SDV4)
	else
	  write(6,801) nat4
	  stop
	endif
        ANG = EL2 * EL1 + EM2 * EM1 + EN2 * EN1
        IF(ANG.LT.0.0) THEN
         EL1 = -EL1
         EM1 = -EM1
         EN1 = -EN1
        ENDIF
        BNORM(1,4) = EL1
        BNORM(2,4) = EM1
        BNORM(3,4) = EN1
C
        EL1 = EL1 + EL2
        EM1 = EM1 + EM2
        EN1 = EN1 + EN2
        CALL NORMAL(EL1,EM1,EN1,EL2,EM2,EN2)
C
      ELSE
C----------------------------------------------------------------------
c      CALL SINGAXS(IPS,N,NAT3,NC8,X,2,XAXIS,YAXIS,EL2,EM2,EN2,IBASE3,S2)
	call findybs(ibapr3,xbase3,ybase3,zbase3,ang1,orig3,'W','C')
	do i=1,3
	  xaxis(2,i)=ybase3(i)
	  yaxis(2,i)=zbase3(i)
	  origin(2,i)= orig3(i)
	enddo
	el2 = xbase3(1)
	em2 = xbase3(2)
	en2 = xbase3(3)
C
      ENDIF
C
      ANG = EL * EL2 +  EM * EM2 +  EN * EN2
      IF(ANG.LT.0.0) THEN
         EL2 = -EL2
         EM2 = -EM2
         EN2 = -EN2
         IF(.NOT.DOUBL2) THEN
c            DO 2038 KL=1,3			! Commented May 23, 2003
c               YAXIS(2,KL) = -YAXIS(2,KL)      ! Commented May 23, 2003
c 2038	    CONTINUE				! Commented May 23, 2003
         ELSE
            DO 2039 KL=1,3
               BNORM(KL,3) = -BNORM(KL,3)
               BNORM(KL,4) = -BNORM(KL,4)
 2039	    CONTINUE
         ENDIF
       ENDIF
C
      IF(.NOT.DOUBL2) THEN
C
c        THETA = 10.0601                        ! bhatta July 2005
c        IF(NAT3.EQ.10) THETA = 12.6921
c         ELS = real(EL2)
c         EMS = real(EM2)
c         ENS = real(EN2)
c         CALL ROTMAT(ELS,EMS,ENS,THETA,R)
c         CALL ROTVEC(R,XAXIS(2,1),XAXIS(2,2),XAXIS(2,3),TEMPX)
c         CALL ROTVEC(R,YAXIS(2,1),YAXIS(2,2),YAXIS(2,3),TEMPY)
c         DO 2040 KJ=1,3			! bhatta July 2005
c            XAXIS(2,KJ) = TEMPX(KJ)
c            YAXIS(2,KJ) = TEMPY(KJ)
c            ORIGIN(2,KJ) = RC8(KJ) - YAXIS(2,KJ) * 4.9122
c 2040    CONTINUE
C
      ELSE
C
      CALL CROSS(YAXIS(2,1),YAXIS(2,2),YAXIS(2,3),EL2,EM2,EN2,
     1XAX(1),XAX(2),XAX(3))
      CALL NORMAL(XAX(1),XAX(2),XAX(3),XAXIS(2,1),XAXIS(2,2),XAXIS(2,3))
c      CALL DIRECX(XAXIS,2,XC1P2,XC82,YAXIS,ybs2,NDX,type3,type4,
c     1cystrns2)    !Sukanya, 12 Nov, 2013
         call direcx(xaxis,2,possblz,el2,em2,en2)
        IF(NDX.EQ.1) THEN
C              TMPBAS = BASE(NOBASE+1)(1:1)
C              BASE(NOBASE+1)(1:1) = BASE(NOBASE+1)(3:3)
C              BASE(NOBASE+1)(3:3) = TMPBAS
           DO 2041 KJ=1,3
              BNORM(KJ,3) = -BNORM(KJ,3)
              BNORM(KJ,4) = -BNORM(KJ,4)
 2041	   CONTINUE
        END IF
       END IF
C
C     FINDING THE X-AXIS OF THE TOP BASE PAIR PLANE
C
       DO 2042  I=1,3
         BM(I) = ORIGIN(2,I) - ORIGIN(1,I)
 2042  CONTINUE
C
c        WRITE(9,*) 'FIRST base-pair NORMAL',EL,EM,EN
        write(9,641) el,em,en
        WRITE(21,'(''BP-ZAXIS(1)'',3F10.6)') EL,EM,EN
641     format('BP-ZAXIS(1)',3f10.6)
        ZAXIS(1,1) = EL
        ZAXIS(1,2) = EM
        ZAXIS(1,3) = EN
        WRITE(9,*) 'SECOND base-pair NORMAL', EL2,EM2,EN2
        WRITE(21,'(''BP-ZAXIS(2)'',3F10.6)')  EL2,EM2,EN2
        ZAXIS(2,1) = EL2
        ZAXIS(2,2) = EM2
        ZAXIS(2,3) = EN2
        LNDEX=LNDEX + 1
	if(ntim.eq.1) then	!  Added bhatta March 30, 2007
          DO 2043 KK1=1,3
            BPNORML(LNDEX,KK1) = ZAXIS(1,KK1)
            MBPNORML(LNDEX,KK1) = (ZAXIS(1,KK1) + ZAXIS(2,KK1)) * 0.5
            BPNORML(LNDEX+1,KK1) = ZAXIS(2,KK1)
            AMFO(LNDEX,KK1) = (ORIGIN(1,KK1) + ORIGIN(2,KK1))*0.5
 2043	  CONTINUE
	endif
      WRITE(9,'(a20,3f10.4)')'MEAN B-PAIR ORIGIN',
     &(AMFO(LNDEX,KK1),KK1=1,3)
C
      WRITE(9,*) '(ORIGIN) CENTER TO CENTER VECTOR _'
      WRITE(9,'(''BM'',3F8.3)')(BM(I),I=1,3)
      WRITE(9,'(''ORIGIN'',3F8.2)')((ORIGIN(J,I),I=1,3),J=1,2)
      IF(DOUBL1) THEN
C
C        WRITE(9,*) ' (1st B.P.) BASE - NORMALS ---'
C 	WRITE(9,'(3F10.6,''  DV ='',F10.7)') (BNORM(M,1),M=1,3),SDV1
C        WRITE(9,'(3F10.6,''  DV ='',F10.7)') (BNORM(M,2),M=1,3),SDV2
C
         CALL BUCK(YAXIS,1,BNORM,PR1,PR2)
         CALL BUCK(XAXIS,1,BNORM,BK1,BK2)
         PROTW(LNDEX) = PR1
C Buckle sign consistent with prop. definition & EMBO convention,
C  10th Aug 1995,  M.B.
C
         BUCAN(LNDEX) = BK1
         CALL BUCK(ZAXIS,1,BPDIR,PR1,PR2)
         OPENAN(LNDEX) = PR2
      ELSE
C       WRITE(9,'(''Base Nornal'',3F10.6,'' DV ='',F10.7)') EL,EM,EN,SD1
      END IF
C
       IF(DOUBL2) THEN
C        WRITE(9,*) ' (2nd. B.P.) BASE - NORMALS---'
C 	WRITE(9,'(3F10.6,''  DV ='',F10.7)') (BNORM(M,3),M=1,3),SDV3
C 	WRITE(9,'(3F10.6,''  DV ='',F10.7)') (BNORM(M,4),M=1,3),SDV4
C
         CALL BUCK(YAXIS,2,BNORM,PR2,PR1)
         CALL BUCK(XAXIS,2,BNORM,BK2,BK1)
         PROTW(LNDEX+1) = PR2
         BUCAN(LNDEX+1) = BK2
         CALL BUCK(ZAXIS,2,BPDIR,PR1,PR2)
         OPENAN(LNDEX+1) = PR2
	ELSE
C     WRITE(9,'(''Base Normal'',3F10.6,'' DV ='',F10.7)')EL2,EM2,EN2,S2
      END IF
C
c Buckle Correction to Base Pair Centers.  ADDED March 29, 2007 from MBU 
c VERSION of NUPARM
c
       CON = 3.14159265/180.0
      IF((BASE(LNDEX)(1:1).EQ."G".AND.BASE(LNDEX+1)(1:1).EQ."A").OR.
     &   (BASE(LNDEX)(1:1).EQ."G".AND.BASE(LNDEX+1)(1:1).EQ."T").OR.
     &   (BASE(LNDEX)(1:1).EQ."C".AND.BASE(LNDEX+1)(1:1).EQ."A").OR.
     &  (BASE(LNDEX)(1:1).EQ."C".AND.BASE(LNDEX+1)(1:1).EQ."T"))THEN
         DZBUCK(LNDEX)  = bproll(LNDEX) - bproll(LNDEX+1)
c	write(*,*) lndex,dzbuck(lndex)
C         DZBUCK(LNDEX)  = BUCAN(LNDEX) - BUCAN(LNDEX+1)   ! With Old Buckle
        ELSE
         DZBUCK(LNDEX)  = bproll(LNDEX+1) - bproll(LNDEX)
c	write(*,*) lndex,dzbuck(lndex)
C         DZBUCK(LNDEX)  = BUCAN(LNDEX+1) - BUCAN(LNDEX)   ! With Old Buckle
        ENDIF
         APROP(LNDEX)  = (PROTW(LNDEX+1) + PROTW(LNDEX))/2.0
         DPROP(LNDEX)  = PROTW(LNDEX+1) - PROTW(LNDEX)
C        DZBUCK(LNDEX)  = BUCAN(LNDEX+1) - BUCAN(LNDEX)
         DOPEN(LNDEX)  = OPENAN(LNDEX+1) - OPENAN(LNDEX)
C           BKAP1 = TAN(BUCAN(LNDEX)*CON*0.5)	! With Old Buckle
C           BKAP2 = TAN(BUCAN(LNDEX+1)*CON*0.5)	! With Old Buckle
           BKAP1 = TAN(bproll(LNDEX)*CON*0.5)
           BKAP2 = TAN(bproll(LNDEX+1)*CON*0.5)
C       WRITE(*,*)BUCAN(LNDEX),BUCAN(LNDEX+1),OPENDS(LNDEX),BKAP1,BKAP2
          DZ1 =  0.25*(OPENDS(LNDEX)) *  BKAP1
          DZ2 =  0.25*(OPENDS(LNDEX+1)) * BKAP2
c	write(*,*) 'BKAP1, BKAP2',bkap1,bkap2,' DZ1, DZ2',dz1,dz2

	write(9,*) 'New Buckle',bproll(lndex),bproll(lndex+1)
	write(9,*) 'New Propeller',bptwst(lndex),bptwst(lndex+1)
	write(9,*) 'New OpenAngle',bptilt(lndex),bptilt(lndex+1)
      WRITE(9,*) 'BP-XAXIS(1)',(XAXIS(1,I),I=1,3)
      WRITE(9,*) 'BP-XAXIS(2)',(XAXIS(2,I),I=1,3)
      WRITE(9,*) 'BP-YAXIS(1)',(YAXIS(1,I),I=1,3)
      WRITE(9,*) 'BP-YAXIS(2)',(YAXIS(2,I),I=1,3)
      WRITE(21,'(''BP-XAXIS(1)'',3F10.6)') (XAXIS(1,I),I=1,3)
      WRITE(21,'(''BP-XAXIS(2)'',3F10.6)') (XAXIS(2,I),I=1,3)
      WRITE(21,'(''BP-YAXIS(1)'',3F10.6)') (YAXIS(1,I),I=1,3)
      WRITE(21,'(''BP-YAXIS(2)'',3F10.6)') (YAXIS(2,I),I=1,3)

	write(9,*) '(ORIGIN) Center to Center Vector'
	do 2052 i=1,3
c	   origin(1,i) = origin(1,i) - (dz1*Zaxis(1,i))
c	   origin(2,i) = origin(2,i) - (dz2*zaxis(2,i))
C   Buckle Correction avoided to make LOCAL to HELICAL transformation exact
C   during MBU visit in June 2015.  bhatta/mb/Prasun
C
	   bm(i) = origin(2,i) - origin(1,i)
	   call normal(bm(1),bm(2),bm(3),bmn(1),bmn(2),bmn(3))
2052 	continue
c	write(*,*) 'BM',bm
	write(9,'(3f9.4)') (bm(i),i=1,3)
C
C  Trial Calculation of basepair OVERLAP as BM.Z1 + BM.Z2
C
      CALL WEGDM(XAXIS,YAXIS,BM,STKINV,NRTCAL,ntim)
C	WRITE(*,*) DXLOC(LNDEX),DYLOC(LNDEX),DZLOC(LNDEX)
      DO 2044 KJ=1,3
        BSCENTR(LNDEX,KJ) = ORIGIN(1,KJ)
        BSCENTR((LNDEX+1),KJ) = ORIGIN(2,KJ)
        HORIG(LNDEX,KJ) = BSCENTR(LNDEX,KJ) - DXLOC(LNDEX) * XAXIS(1,KJ)
     1                  - DYLOC(LNDEX) * YAXIS(1,KJ)
 2044 CONTINUE
      IF(STKINV) THEN 
  	 CALL WEGDM(XAXIS,YAXIS,BM,STKINV,NRTCAL,ntim)
      END IF
      CALL GLOBAL(ORIGIN,XAXIS,YAXIS,STKINV)
C ----------------------------------------------------------------------
C      CLEARING THE ARRAYS CONTAINING COORDINATES TO 0.0 BEFORE GOING
C              TO THE NEXT base-pair STEP.
C ----------------------------------------------------------------------
      DO 2045 K=1,3
          XC8(K) = 0.0
          XC6(K) = 0.0
          RC8(K) = 0.0
          RC6(K) = 0.0
          PROJ(K) = 0.0
          ORIG(K) = 0.0
          XYZ1(K) = 0.0
          XYZ2(K) = 0.0
          XYZ3(K) = 0.0
          XYZ4(K) = 0.0
          SM(K) = 0.0D0
          BM(K) = 0.0D0
          BMN(K) = 0.0D0
 2045 CONTINUE
      DO 2046 K=1,15
           NSP(K) = 0
           NSP1(K) = 0
           N(K) = 0
           N1(K) = 0
           DVA(K) = 0.0
           DO 2047 J=1,3
              X(J,K) = 0.0
              X1(J,K) = 0.0
 2047	   CONTINUE
 2046	CONTINUE
        WRITE(9,'(3(''****************************''))')
      RETURN
      END
C ********************************************************************
 
      SUBROUTINE CROSS(A,B,C,A1,B1,C1,ALPHA,BETA,GAMA)
      DOUBLE PRECISION A,B,C,ALPHA,BETA,GAMA, A1,B1,C1
      ALPHA=B*C1-C*B1
      BETA=C*A1-A*C1
      GAMA=A*B1-B*A1
      AMOD=DSQRT(ALPHA*ALPHA+BETA*BETA+GAMA*GAMA)
      ALPHA=ALPHA/AMOD
      BETA=BETA/AMOD
      GAMA=GAMA/AMOD
      RETURN
      END
 
C --------------------------------------------------------------------C
C  NORMALISING SUBROUTINE
C --------------------------------------------------------------------
      SUBROUTINE NORMAL(A,B,C,X,Y,Z)
      DOUBLE PRECISION A,B,C, X,Y,Z, ALEN,ALENGT
C 
      ALEN =  A * A + B * B + C * C
C 
      ALENGT = DSQRT(ALEN)
C 
      IF (ABS(ALENGT) .LT. 1.0E-6) THEN
         X=0.D0
         Y=0.D0
         Z=0.D0
      ELSE
         X = A/ALENGT
         Y = B/ALENGT
         Z = C/ALENGT
      ENDIF
C 
      RETURN
      END
C --------------------------------------------------------------------
C    SUBROUTINE FOR TIP AND INCL. CALC. W.R.T. LOCAL HELIX AXIS &
C     WTILT AND WROLL CALC. W.R.T. THE MEAN DOUBLET AXES.
C --------------------------------------------------------------------
      SUBROUTINE WEGDM(XAXIS,YAXIS,BMN,STKINV,NRTCAL,ntim)
	parameter ( nrs = 9999 )
	include 'parameters.h'
      LOGICAL STKINV,RETCAL
      DIMENSION XAXIS(2,3),YAXIS(2,3),XM(3),YM(3),ZAXIS(3),BMN(3)
      COMMON /BPZAXS/BPNORML(nrs,3),BPNORMLX(nrs,3),BPNORMLY(nrs,3),
     1     MBPNORML(nrs,3),AMFO(nrs,3)
      COMMON /REPLY/ANSOR,ANSC1,ANSBPN,ANSDB,ANSrnt,ANSHO,ansnrm,answw,
     1    anspp,anstor,ansovl,anshlx
      COMMON /ATOMS/IBASE1(10),IBASE2(10),IBASE3(10),IBASE4(10),NB,
     1 ibapr(nrs)
      DOUBLE PRECISION XAXIS,YAXIS,ZAXIS,XM,YM,ZM1,ZM2,ZM3,Y1,Y2,Y3,
     1   Y12,Y22,Y32,X11,X12,X13,BMN,TL2,RL2,WTILT,WROLL,Y11,Y13,
     2   BPNORML,BPNORMLX,BPNORMLY,MBPNORML,amfo
      character*1 ansdb,ansbpn,ansnrm,ansor,ansc1,ansrnt,ansho,answw,
     1 anssng,anszp,anstrj,anspp,anstor,angl,ansovl,anshlx
C
      CONV=180.0/3.141592654
C
      XM(1) = XAXIS(1,1) - XAXIS(2,1)
      XM(2) = XAXIS(1,2) - XAXIS(2,2)
      XM(3) = XAXIS(1,3) - XAXIS(2,3)
      YM(1) = YAXIS(1,1) - YAXIS(2,1)
      YM(2) = YAXIS(1,2) - YAXIS(2,2)
      YM(3) = YAXIS(1,3) - YAXIS(2,3)
      CALL CROSS(XM(1),XM(2),XM(3),YM(1),YM(2),YM(3),ZM1,ZM2,ZM3)
      CALL NORMAL(ZM1,ZM2,ZM3,ZAXIS(1),ZAXIS(2),ZAXIS(3))
      DO 2000 KJ=1,3
          ZAXISH(KJ,LNDEX) = ZAXIS(KJ)
 2000 CONTINUE
      CALL CROSS(XAXIS(1,1),XAXIS(1,2),XAXIS(1,3),ZAXIS(1),ZAXIS(2),ZAXI
     1S(3),Y1,Y2,Y3)
      CALL CROSS(XAXIS(2,1),XAXIS(2,2),XAXIS(2,3),ZAXIS(1),ZAXIS(2),ZAXI
     1S(3),Y12,Y22,Y32)
      ANG = (Y1 * Y12 + Y2 * Y22 + Y3 * Y32)/((Y1*Y1+Y2*Y2+Y3*Y3)*
     1   (Y12*Y12+Y22*Y22+Y32*Y32))
      ANGLE = ACOS(ANG) * CONV

c      IF(ABS(ANGLE).GT.90.0.AND.(.NOT.STKINV)) STKINV = .TRUE.
      IF(ABS(ANGLE).GT.120.0.AND.(.NOT.STKINV)) STKINV = .TRUE.

      CALL CROSS(Y1,Y2,Y3,Y12,Y22,Y32,ZM1,ZM2,ZM3)
      ANG = ZM1 * ZAXIS(1) + ZM2 * ZAXIS(2) + ZM3 * ZAXIS(3)
C
      IF(ANG.LT.0) ANGLE = -ANGLE
C      IF(NRTCAL.EQ.0) 
 	TWISTH(LNDEX) = ANGLE
        WRITE(9,99)
 99   FORMAT(/)
        WRITE(9,*) 'Z*-AXIS',(ZAXIS(I),I=1,3)
        WRITE(21,'(''Z*-Axis or Local Helix Axis'',3F10.6)') 
     1 (ZAXIS(I),I=1,3)
      WRITE(21,'(''Local Helix Origin'',3f10.3)')
     &  (AMFO(LNDEX,KK1),KK1=1,3)
      WRITE(21,*) '****** Step no. ******'
C      END IF
      AINC = 0.0
      TIP1 = 0.0
      TSLIDE = 0.0
      DO 2001 I=1,3
         AINC = AINC + YAXIS(1,I) * ZAXIS(I)
         TIP1 = TIP1 + XAXIS(1,I) * ZAXIS(I)
         TSLIDE = TSLIDE + BMN(I) * ZAXIS(I)
 2001 CONTINUE
C
      DZL(LNDEX) = TSLIDE
      AINCL =    ASIN(AINC) * CONV
      TIPL =  -ASIN(TIP1) * CONV
      DIETLT(LNDEX) = AINCL
      DIEROL(LNDEX) = TIPL

C
      XM(1) = XAXIS(1,1) + XAXIS(2,1)
      XM(2) = XAXIS(1,2) + XAXIS(2,2)
      XM(3) = XAXIS(1,3) + XAXIS(2,3)
      YM(1) = YAXIS(1,1) + YAXIS(2,1)
      YM(2) = YAXIS(1,2) + YAXIS(2,2)
      YM(3) = YAXIS(1,3) + YAXIS(2,3)
      CALL CROSS(XM(1),XM(2),XM(3),YM(1),YM(2),YM(3),ZM1,ZM2,ZM3)
      CALL NORMAL(ZM1,ZM2,ZM3,ZAXIS(1),ZAXIS(2),ZAXIS(3))
      CALL NORMAL(XM(1),XM(2),XM(3),X11,X12,X13)
      CALL NORMAL(YM(1),YM(2),YM(3),Y11,Y12,Y13)
C
      XM(1) = X11
      XM(2) = X12
      XM(3) = X13
      YM(1) = Y11
      YM(2) = Y12
      YM(3) = Y13

	if(ntim.eq.1) then	! Added bhatta March 30, 2007
	  do 2200 kj=1,3
	     bpnorml(lndex,kj) = zaxis(kj)
	     bpnormlx(lndex,kj) = xm(kj)
	     bpnormly(lndex,kj) = ym(kj)
 2200	  continue
	endif
      angxmym=acos(x11*y11+x12*y12+x13*y13)*conv
      write(9,*) 'Angle between Xm & Ym=',angxmym
      TL2 = 0.0
      RL2 = 0.0
      CNEW = 0.0
      TSLIDE = 0.0
      SLIDE = 0.0
      WRITE(9,*) 'MEAN Z-AXIS ',(ZAXIS(I),I=1,3)
      DO 2002 I=1,3
         TL2 = TL2 - YAXIS(1,I) * ZAXIS(I)
         RL2 = RL2 + XAXIS(1,I) * ZAXIS(I)
         CNEW = CNEW + XM(I) * BMN(I)
         TSLIDE = TSLIDE + BMN(I) * ZAXIS(I)
         SLIDE = SLIDE + BMN(I) * YM(I)
 2002 CONTINUE
C
      WTILT =  2.0 * DASIN(TL2) * CONV
      WROLL =  2.0 * DASIN(RL2) * CONV
C
      TILTL(LNDEX) = WTILT
      ROLLL(LNDEX) = WROLL
      SLY(LNDEX) = SLIDE
      SLX(LNDEX) = CNEW
      SLZ(LNDEX) = TSLIDE
	write(9,*) 'Tilt, roll etc. just after calcn',wtilt,wroll,
     1 slide,cnew,tslide
      CALL CROSS(XAXIS(1,1),XAXIS(1,2),XAXIS(1,3),ZAXIS(1),ZAXIS(2),ZAXI
     1S(3),Y1,Y2,Y3)
      CALL CROSS(XAXIS(2,1),XAXIS(2,2),XAXIS(2,3),ZAXIS(1),ZAXIS(2),ZAXI
     1S(3),Y12,Y22,Y32)
      ANG = Y1 * Y12 + Y2 * Y22 + Y3 * Y32
      ANGLE = ACOS(ANG) * CONV
C
      CALL CROSS(Y1,Y2,Y3,Y12,Y22,Y32,ZM1,ZM2,ZM3)
      ANG = ZM1 * ZAXIS(1) + ZM2 * ZAXIS(2) + ZM3 * ZAXIS(3)
      IF(ANG.LT.0.0) ANGLE = - ANGLE
C
C      IF(NRTCAL.EQ.0)  
	TWISTL(LNDEX) = ANGLE
      CALL SLIDCAL
      IF(STKINV) THEN
  	  DO 3000 KK=1,3
	    XAXIS(2,KK) = -XAXIS(2,KK)
	    YAXIS(2,KK) = -YAXIS(2,KK)
 3000	  CONTINUE
      END IF
C
      IF(NRTCAL.NE.0) STKINV = .FALSE.
      NRTCAL = NRTCAL + 1
C
      RETURN
      END
C ----------------------------------------------------------------------
C       CALCULATION OF PROPELLER TWIST AND BUCKLE FROM THE BASE-NORMALS.
C ----------------------------------------------------------------------
      SUBROUTINE BUCK(AXIS,I,BNORML,VALUE2,VALUE)
      DOUBLE PRECISION AXIS,PROP,BNORML,P1,P2,P4,P5
      DIMENSION AXIS(2,3),PROP(3,2),BNORML(3,4),P1(3),P2(3),P4(3),P5(3)
C
      CONV=180.0/3.141592654
      J=1
      IF(I.EQ.2) J=3
      J1 = J + 1
C
C  CALCULATION OF PARAMETERS TAKING THE PROJECTION BY METHOD OF
C  "GRAM-SMGITH".
C --------------------------------------------------------------------
       PROJEC = 0.0
       PROJE2 = 0.0
       PRTW = 0.0
       SIGN = 0.0
       DO 2000 KJ=1,3
          PROJEC = PROJEC + AXIS(I,KJ) * BNORML(KJ,J)
          PROJE2 = PROJE2 + AXIS(I,KJ) * BNORML(KJ,J1)
 2000  CONTINUE
       DO 2001 KJ=1,3
          PROP(KJ,1) = BNORML(KJ,J) - PROJEC * AXIS(I,KJ)
          PROP(KJ,2) = BNORML(KJ,J1) - PROJE2 * AXIS(I,KJ)
 2001  CONTINUE
       CALL NORMAL(PROP(1,1),PROP(2,1),PROP(3,1),P1(1),P1(2),P1(3))
       CALL NORMAL(PROP(1,2),PROP(2,2),PROP(3,2),P2(1),P2(2),P2(3))
       DO 2002 KJ=1,3
          PRTW = PRTW + P1(KJ) * P2(KJ)
 2002  CONTINUE
       VALUE = ACOSF(PRTW)
C --------------------------------------------------------------------
C   CALCULATION OF PARAMETERS BY TAKING CROSS PRODUCT.
C---------------------------------------------------------------------
      CALL CROSS(AXIS(I,1),AXIS(I,2),AXIS(I,3),BNORML(1,J),BNORML(2,J),
     1BNORML(3,J),P4(1),P4(2),P4(3))
      CALL CROSS(AXIS(I,1),AXIS(I,2),AXIS(I,3),BNORML(1,J1),BNORML(2,J1)
     1,BNORML(3,J1),P5(1),P5(2),P5(3))
      PRTW = 0.0
      DO 2003 KJ=1,3
         PRTW = PRTW + P4(KJ) * P5(KJ)
 2003 CONTINUE
      VALUE2 = ACOSF(PRTW)
      SIGN = 0.0
      CALL CROSS(P4(1),P4(2),P4(3),P5(1),P5(2),P5(3),P1(1),P1(2),P1(3))
      DO 2004 KJ=1,3
         SIGN = SIGN + P1(KJ) * AXIS(I,KJ)
 2004 CONTINUE
      IF(SIGN.GT.0.0) VALUE2 = - VALUE2
      RETURN
      END
C----------------------------------------------------------------------
C  SUBROUTINE TO FIND OUT DIRECTION TOWARD MAJOR GROOVE, TO FIX X-AXIS.
C----------------------------------------------------------------------
      SUBROUTINE DIRECX2(AXIS,I,XC8,XN9,Y,YBS,IND,typ1,typ2,bptyp)
      DOUBLE PRECISION AXIS,VEC,VECN,Y,YBS
      DIMENSION AXIS(2,3),XC8(3),XN9(3),VEC(3),VECN(3),Y(2,3),YBS(2,3)
      character*1 typ1,typ2,bptyp

      IND = 0
      DO 2000 JK = 1,3
         VEC(JK) = XC8(JK) - XN9(JK)
 2000 CONTINUE
      CALL NORMAL(VEC(1),VEC(2),VEC(3),VECN(1),VECN(2),VECN(3))
      ANG = 0.0
      DO 2001 JK = 1,3
         IF((typ1.EQ.'S'.or.typ1.eq.'s'.or.typ1.eq.'z')   !Sukanya, 12 Nov, 2013
     1.and.(bptyp.eq.'T'))then
           YBS(1,JK)=-YBS(1,JK)
         END IF 
            ANG = ANG + AXIS(I,JK) * YBS(1,JK)
 2001 CONTINUE
      IF(ANG.LT.0.0) THEN
        IND = 1
           DO 2002 KJ=1,3
                      AXIS(I,KJ) = -AXIS(I,KJ)
                      Y(I,KJ) = -Y(I,KJ)
 2002      CONTINUE
      ENDIF
      RETURN
      END
C ---------------------------------------------------------------------
C
C     SUBROUTINE FOR CALCULATION OF HELICAL (GLOBAL) PARAMETERS.
C
C ---------------------------------------------------------------------
        SUBROUTINE GLOBAL(ORIG,XAXIS,YAXIS,STKINV)
	parameter ( nrs = 9999 )
	include 'parameters.h'
        DOUBLE PRECISION XAXIS,YAXIS
	 LOGICAL STKINV
        DIMENSION ORIG(2,3),XAXIS(2,3),YAXIS(2,3)
C
        CONV = 180.0/3.1415926
C
        HREPIT = ORIG(1,3) - ORIG(2,3)
        IF(HREPIT.LT.0.0) HREPIT = -HREPIT
        HG(LNDEX) = HREPIT
C
        CALL DXDY(ORIG,YAXIS,1,DX1,DY1)
        CALL DXDY(ORIG,YAXIS,2,DX2,DY2)
        DXG(LNDEX) = DX1
        DYG(LNDEX) = DY1
        DXG(LNDEX+1) = DX2
        DYG(LNDEX+1) = DY2
C
        TILT = ASIN(YAXIS(1,3)) * CONV
        ROLL = ASIN(-XAXIS(1,3)) * CONV
        TILT2 = ASIN(YAXIS(2,3)) * CONV
        ROLL2 = ASIN(-XAXIS(2,3)) * CONV
        TILTG(LNDEX) = TILT
        ROLLG(LNDEX) = ROLL
        TILTG(LNDEX+1) = TILT2
        ROLLG(LNDEX+1) = ROLL2
C
        DIM = YAXIS(1,1) *YAXIS(1,1) + YAXIS(1,2) * YAXIS(1,2)
        DIM = SQRT(DIM)
        Y1Z1 = YAXIS(1,1) / DIM
        Y1Z2 = YAXIS(1,2) / DIM
C
        DIM = YAXIS(2,1) * YAXIS(2,1) + YAXIS(2,2) * YAXIS(2,2)
        DIM = SQRT(DIM)
        Y2Z1 = YAXIS(2,1) / DIM
        Y2Z2 = YAXIS(2,2) / DIM
        ANG = Y1Z1 * Y2Z1 + Y1Z2 * Y2Z2
        ANGLE = ACOS(ANG) * CONV
        TWISTG(LNDEX) = ANGLE
        ZCOM = Y1Z1 * Y2Z2 - Y1Z2 * Y2Z1
	 IF(ZCOM.LT.0.0) TWISTG(LNDEX) = -TWISTG(LNDEX)
        RETURN
        END
C
        SUBROUTINE DXDY(ORIGIN,AXIS,I,DX,DY)
        DOUBLE PRECISION AXIS
C       REAL LY2MY2,NUM
        DIMENSION ORIGIN(2,3),AXIS(2,3)
C
        Y2LY2M = AXIS(I,1) * AXIS(I,1) + AXIS(I,2) * AXIS(I,2)
        GRAD = AXIS(I,2) / AXIS(I,1)
        ANUM = ORIGIN(I,2) - GRAD * ORIGIN(I,1)
C
C	write(6,*) 'y2ly2m =',y2ly2m,' and axis',(axis(i,k),k=1,3)
        DX = -AXIS(I,1) * ANUM /(SQRT(Y2LY2M))
C
        X1 = -ANUM * AXIS(I,1) * AXIS(I,2) / Y2LY2M
        Y1 = AXIS(I,1) * AXIS(I,1) * ANUM / Y2LY2M
C
        DY = (X1-ORIGIN(I,1)) * (X1-ORIGIN(I,1)) + (Y1-ORIGIN(I,2)) *
     1 (Y1-ORIGIN(I,2))
C     	write(6,*) 'dy before sqrt is',dy
        DY = SQRT(DY)
        EL = (X1 - ORIGIN(I,1))/DY
        EM = (Y1 - ORIGIN(I,2))/DY
        ANG = EL * AXIS(I,1) + EM * AXIS(I,2)
        IF(ANG.GT.0.0) DY = -DY
        RETURN
        END
C
C-----------------------------------------------------------------------
      SUBROUTINE WRITING(FILENM,TITLE1,INPUT,NTIM,NRUN,ANORIN,filprt,
     1 hlxfl)
	parameter ( nrs = 9999 )
	include 'parameters.h'
      DOUBLE PRECISION XAXIS,YAXIS,ZAXIS,XM,YM,ZM1,ZM2,ZM3,Y1,Y2,Y3,
     1   Y12,Y22,Y32,X11,X12,X13,BMN,TL2,RL2,WTILT,WROLL,Y11,Y13,amfo
        COMMON /BASINF/BASE(nrs),pair(3,nrs),prtype(nrs)
      DOUBLE PRECISION BPNORML,BPNORMLX,BPNORMLY,MBPNORML
      COMMON /BPZAXS/BPNORML(nrs,3),BPNORMLX(nrs,3),BPNORMLY(nrs,3),
     1     MBPNORML(nrs,3),AMFO(nrs,3)
      COMMON /GLOREO/DCRS1,DCRS2,DCRS3,ANGLE,DX1,DY1,INFBP(2,nrs)
      COMMON /REPLY/ANSOR,ANSC1,ANSBPN,ANSDB,ANSrnt,ANSHO,ansnrm,answw,
     1    anspp,anstor,ansovl,anshlx
      DIMENSION ANGST(nrs,nrs),TITLE1(20),ihlx(2,nrs),hlxpr(nrs),
     1 hlxor(nrs)
	common /nameorg/secnmo(nrs)
      CHARACTER *80 TITLE,TITLE1,filprt*20
      character*40 filenm,input,hlxfl
      CHARACTER *4 ANORIN
      CHARACTER*3 BASE,bpadd,hlxpr,secnmo
      CHARACTER*1 ANSOR,ANSC1,ANSrnt,ANSHO,ANSDB,ANSBPN,pair,ansnrm,
     1            prtype,answw,anspp,anstor,ansovl,hlxor,anshlx
C
      CONV=180.0/3.141592654
      zero = 0.0
      one = 1.0
      if((anshlx.eq.'Y').or.(anshlx.eq.'y'))then
        open(unit=27,file=hlxfl)
        do IH=1,9999
          read(27,178,end=179)ihlx(1,ih),ihlx(2,ih),hlxpr(ih),hlxor(ih)
          write(*,178)ihlx(1,ih),ihlx(2,ih),hlxpr(ih),hlxor(ih)
        enddo
179     continue
      endif
      infopr=ih
178   format(I5,1X,I5,1X,A3,1X,A1)
C 

      IF(NTIM .EQ. 1)THEN
C 
        KKKK=INDEX(FILENM(1:40),':')
C 
        IF (KKKK .EQ. 0) THEN
         KKKK=1
        ELSE
         KKKK=KKKK+1
        ENDIF
C 
        KKKKK=INDEX(FILENM(KKKK:40),'.')
        KKKKK=KKKKK+KKKK
C 
        TITLE = FILENM(KKKK:KKKKK-1)//'prm'
        filprt = filenm(kkkk:kkkkk-1)
C
      OPEN(UNIT=11,FILE=TITLE)
      WRITE(11,300) FILENM
c      WRITE(6,400) TITLE
      WRITE(8,400) TITLE
      title = filenm(kkkk:kkkkk-1)//'out'
      open(unit=52,file=title)
 300  FORMAT(1X,'COORDINATES INPUT TAKEN FROM ::',A40)
 400  FORMAT(/' Results of the calculation are written into: '/1X,A40)
 401  format(a80)
c       DO 2000 II=1,KKTI
c          WRITE(11,'(A80)') TITLE1(II)
c 2000 	CONTINUE
c	do ii=1,9999
c	  read(52,401,end=998) title
c	  write(11,401) title
c	enddo
 998	continue
        IF(ANSBPN.NE.'Y') WRITE(11,61)INPUT
 61   FORMAT(/1X,'RESIDUE SERIAL NUMBERS TAKEN FROM ::',A40)
	write(11,997)
 997	format(' Basepairing information has been calculated using', 
     1 'BPFIND software, Das et al. 2006, J. Biomol. Struct. Dynam. 24,
     2 149'/'Parameters are defined in (i) Bansal et al (1995) CABIOS',
     3 '11, 281'/25x,'(ii) Mukherjee et al (2006) J. Comp. Aided Mol.'
     4 ,' Des. 20, 629')
 
	if(ansnrm.eq.'Y') write(11,51)
 51	format(/2x,'* Base normals are calculated by successive cross-',
     1     'products *')
	if(ansnrm.ne.'Y') write(11,52)
 52	format(/2x,'* Base normals are calculated by least-square fit *')
        IF(ANSOR.EQ.'Y') WRITE(11,71)
 71   FORMAT(/2X,'*Basepair origin is at the center of gravity*')
 
        IF(ANSC1.EQ.'Y') THEN
          WRITE(11,81)
 81   FORMAT(/2X,'The basepair Y-axis is along C1''--C1'' direction.'/)
        ELSE
          WRITE(11,82)
 82   FORMAT(/2X,'The basepair Y-axis is along C6--C8 direction.'/)
        END IF
      END IF
C-----------------------------------------------------------------------
      IF (NTIM .EQ. 2)WRITE(11,221)
 221  FORMAT(//2X,'Single strand parameters for Strand 1 of duplex:')
      IF (NTIM .EQ. 3)WRITE(11,222)
 222  FORMAT(//2X,'Single strand parameters for Strand 2 of duplex:')
 193  FORMAT(2x,'Single Strand Parameters with New Definition of',
     1 'Base fixed axes:')
      if(ansovl.eq.'Y') then
        write(11,4101)
      else
        WRITE(11,1)
      endif
      WRITE(11,199)
      DO K=1,LNDEX+1
	 if(iachar(pair(1,k)).eq.0) pair(1,k) = 'W'
	 if(iachar(pair(2,k)).eq.0) pair(2,k) = 'W'
	 if(iachar(pair(3,k)).eq.0) pair(3,k) = 'C'
      enddo
      DO 2001 K=1,LNDEX
	if(bproll(k).ne.0.and.bproll(k+1).ne.0) then 
	  dzbuck(k)=bproll(k+1) - bproll(k)
	endif
        if(secnmo(infbp(1,k)).eq.'THY') base(k)(1:1)='T'
        if(secnmo(infbp(2,k)).eq.'THY') base(k)(3:3)='T'
       if(ansovl.eq.'Y') then
        WRITE(11,1410) INFBP(1,K),BASE(K),INFBP(2,K), TILTL(K),ROLLL(K),
     1   TWISTL(K),SLX(K),SLY(K),SLZ(K),dzbuck(k),ovrlp(k),pair(1,k),
     2   pair(2,k),pair(3,k),filprt
         else
        if(infbp(1,k+1).eq.0) then
           tiltl(k)=zero/zero
           rolll(k)=zero/zero
           twistl(k)=zero/zero
           slx(k) = zero/zero
           sly(k)=zero/zero
           slz(k) = zero/zero
           dzbuck(k) = zero/zero
         endif
         if((anshlx.eq.'Y').or.(anshlx.eq.'y'))then
           do ih=1,infopr
             if((ihlx(1,ih).ne.0).and.(infbp(1,k).eq.ihlx(1,ih)))then
         WRITE(11,410) INFBP(1,K),BASE(K),INFBP(2,K), TILTL(K),ROLLL(K),
     1   TWISTL(K),SLX(K),SLY(K),SLZ(K),dzbuck(k),pair(1,k),pair(2,k),
     2   pair(3,k),filprt
             endif
           enddo
         else
         WRITE(11,410) INFBP(1,K),BASE(K),INFBP(2,K), TILTL(K),ROLLL(K),
     1   TWISTL(K),SLX(K),SLY(K),SLZ(K),dzbuck(k),pair(1,k),pair(2,k),
     2   pair(3,k),filprt
         endif

       endif
 2001 CONTINUE
 410 	format('LC',1X,I4,2X,A3,2X,I4,7(1X,F7.2),1x,a1,':',a1,2x,a1,2x,
     1         a20)
1410 	format('LC',1X,I4,2X,A3,2X,I4,8(1X,F7.2),1x,a1,':',a1,2x,a1,2x,
     1         a20)
 420    format('HL',1X,I4,2X,A3,2X,I4,6(1X,F7.2),1x,a1,':',
     1  a1,2x,a1,2x,a20)
 423 	format('HL',1X,I4,2X,A3,2X,I4,7(1X,F7.2),1x,a1,':',
     1  a1,2x,a1,2x,a20)
 430 	format('BP',1X,I4,2X,A3,2X,I4,7(1X,F7.2),1x,a20)
 440 	format('GL',1X,I4,2X,A3,2X,I4,6(1X,F7.2),1x,a20)
 431	format('BL',1x,i4,2x,a3,2x,i4,6(1x,f7.2),1x,a1,':',a1,2x,a1,2x,
     1         a20)
 434	format('TL',1x,i4,2x,a3,2x,i4,6(1x,f7.2),1x,a1,':',a1,2x,a1,2x,
     1         a20)
 435	format('FL',1x,i4,2x,a3,2x,i4,6(1x,f7.2),1x,a1,':',a1,2x,a1,2x,
     1         a20)
 432	format('TP',1x,i4,2x,a3,2x,i4,6(1x,f7.2),1x,a1,':',a1,2x,a1,2x,
     1         a20)
 433	format('BF',1x,i4,2x,a3,2x,i4,6(1x,f7.2),1x,a1,':',a1,2x,a1,2x,
     1         a20)
 441  FORMAT('GL',1X,I4,2X,A3,2X,I4,2F8.2,8X,2F8.2,9x,a20)
      CALL STATIS(TILTL,LNDEX,AVTL,STDTL)
      CALL STATIS(ROLLL,LNDEX,AVRL,STDRL)
      CALL STATIS(TWISTL,LNDEX,AVTW,STDTW)
      CALL STATIS(SLX,LNDEX,AVSX,STDSX)
      CALL STATIS(SLY,LNDEX,AVSY,STDSY)
      CALL STATIS(SLZ,LNDEX,AVSZ,STDSZ)
      call statis(dzbuck,lndex,avcp,stdcp)
      if(ansovl.eq.'Y') then
        call statis(ovrlp,lndex,avovl,stdovl)
      else
        avovl=0.0
        stdovl=0.0
      endif
      CALL STATIS(PROTW,(LNDEX+1),AVPR,STDPR)
      CALL STATIS(BUCAN,(LNDEX+1),AVBU,STDBU)
      CALL STATIS(OPENAN,(LNDEX+1),AVAN,STDAN)
      CALL STATIS(OPENDS,(LNDEX+1),AVDS,STDDS)
      CALL STATIS(OPENC1,(LNDEX+1),AVC1,STDC1)
      CALL STATIS(ANGBN9,(LNDEX+1),AVBN9,STDBN9)
      CALL STATIS(ANGBN1,(LNDEX+1),AVBN1,STDBN1)
      WRITE(11,410) INFBP(1,(LNDEX+1)),BASE(LNDEX+1),INFBP(2,(LNDEX+1))
       WRITE(11,199)
       if(ansovl.eq.'Y') then
       WRITE(11,110) AVTL,AVRL,AVTW,AVSX,AVSY,AVSZ,avcp,avovl
       WRITE(11,111) STDTL,STDRL,STDTW,STDSX,STDSY,STDSZ,stdcp,stdovl
       else
       WRITE(11,110) AVTL,AVRL,AVTW,AVSX,AVSY,AVSZ,avcp
       WRITE(11,111) STDTL,STDRL,STDTW,STDSX,STDSY,STDSZ,stdcp
       endif 
       WRITE(11,199)
      if(ntim.eq.2.or.ntim.eq.3) then
        WRITE(11,193)
        WRITE(11,1)
        WRITE(11,199)
        do 3001 k=1,lndex
	  write(11,410) infbp(1,k),base(k),infbp(2,k),tilts(k),rolls(k),
     1        twists(k),slxs(k),slys(k),dzls(k)
 3001	continue
      WRITE(11,410) INFBP(1,(LNDEX+1)),BASE(LNDEX+1),INFBP(2,(LNDEX+1))
       WRITE(11,199)
      CALL STATIS(TILTs,LNDEX,AVTL,STDTL)
      CALL STATIS(rolls,LNDEX,AVrl,STDrl)
      CALL STATIS(twists,LNDEX,AVtw,STDtw)
      CALL STATIS(slxs,LNDEX,AVsx,STDsx)
      CALL STATIS(slys,LNDEX,AVsy,STDsy)
      CALL STATIS(dzls,LNDEX,AVsz,STDsz)
       WRITE(11,110) AVTL,AVRL,AVTW,AVSX,AVSY,AVSZ
       WRITE(11,111) STDTL,STDRL,STDTW,STDSX,STDSY,STDSZ
       WRITE(11,199)
      endif
 199  FORMAT(' ----------------------------------------------------------
     1-----------')
 999   FORMAT(//)
 110   FORMAT(11X,'Av.',4X,11F8.2)
 111   FORMAT(9X,'Std.Dev. ',10F8.2)
73     FORMAT(11X,'Av.',4X,11F8.2)
74     FORMAT(9X,'Std.Dev. ',10F8.2)
  1   FORMAT(2X,'Local Step Parameters: '/22X,'Tilt    Roll   Twist    S
     1hift   Slide   Rise    Cup')
4101   FORMAT(2X,'Local Step Parameters: '/22X,'Tilt    Roll   Twist    
     1Shift   Slide   Rise    Cup  Overlap')
  3   FORMAT(/,2X,'Local Helical Parameters:',/20X,'Inclin.   Tip    Twi
     1st',5X,'dx',6X,'dy',6X,'dz')
 100  FORMAT(1X,I2,2X,A3,2X,I2,11(1X,F7.2))
 101  FORMAT(1X,I2,2X,A3,2X,I2,2F8.2,8X,2F8.2)
        WRITE(11,3)
      WRITE(11,199)
      DO 2002 K=1,LNDEX
        if(infbp(1,k+1).eq.0) then
           dietlt(k)=zero/zero
           dierol(k)=zero/zero
           twisth(k)=zero/zero
           dxloc(k) = zero/zero
           dyloc(k)=zero/zero
           dzloc(k) = zero/zero
         endif
        if((anshlx.eq.'Y').or.(anshlx.eq.'y'))then
           do ih=1,infopr
             if((ihlx(1,ih).ne.0).and.(infbp(1,k).eq.ihlx(1,ih)))then
          WRITE(11,420) INFBP(1,K),BASE(K),INFBP(2,K),DIETLT(K),
     1    DIEROL(K),TWISTH(K),DXLOC(K),DYLOC(K),DZLOC(K),
     2    pair(1,k),pair(2,k), pair(3,k),filprt
             endif
           enddo
         else
          WRITE(11,420) INFBP(1,K),BASE(K),INFBP(2,K),DIETLT(K),
     1    DIEROL(K),TWISTH(K),DXLOC(K),DYLOC(K),DZLOC(K),
     2    pair(1,k),pair(2,k), pair(3,k),filprt
         endif

 2002 CONTINUE
      CALL STATIS(DIETLT,LNDEX,AVTL,STDTL)
      CALL STATIS(DIEROL,LNDEX,AVRL,STDRL)
      CALL STATIS(TWISTH,LNDEX,AVTW,STDTW)
      CALL STATIS(DXLOC,LNDEX,AVSX,STDSX)
      CALL STATIS(DYLOC,LNDEX,AVSY,STDSY)
      CALL STATIS(DZLOC,LNDEX,AVSZ,STDSZ)
      WRITE(11,420) INFBP(1,(LNDEX+1)),BASE(LNDEX+1),INFBP(2,(LNDEX+1))
      WRITE(11,199)
        WRITE(11,110) AVTL,AVRL,AVTW,AVSX,AVSY,AVSZ
        WRITE(11,111) STDTL,STDRL,STDTW,STDSX,STDSY,STDSZ
      WRITE(11,199)
C-----------------------------------------------------------------------
       IF(ANSDB.EQ.'Y'.AND.NTIM.LE.1) THEN
        WRITE(11,5)
  5     FORMAT(/2X,'Intra base-pair Parameters: '/21X,'Prop.  Buckle  Op
     1enan*    Glycosidic    C1..C1  C8..C6   ')
        WRITE(11,199)
         DO 2004 K=1,(LNDEX+1)
       if(INFBP(2,K).ne.0)then
       WRITE(11,430)INFBP(1,K),BASE(K),INFBP(2,K),PROTW(K),BUCAN(K),
     1 OPENAN(K),ANGBN9(K),ANGBN1(K),OPENC1(K),OPENDS(K),filprt
       endif
 2004	CONTINUE
         WRITE(11,199)
         WRITE(11,110) AVPR,AVBU,AVAN,AVBN9,AVBN1,AVC1,AVDS
         WRITE(11,111) STDPR,STDBU,STDAN,STDBN9,STDBN1,STDC1,STDDS
         WRITE(11,199)
	 WRITE(11,227)
 227  FORMAT(' *OPENAN is the angle between C6 (or C8)..C1'' vectors 
     1in the basepair,',/,'  projected on the X-Y plane. It has a value
     2 of ~20 in a W-C basepair.')
         write(11,2016)
2016    format(/2x,'Intra Base-pair Parameters following Local Parameter
     1 Definition'/21x,'Buckle  Open    Propel   Stagger Shear Stretch')
     	write(11,199)
         DO 2015 K=1,(LNDEX+1)
	   if(iachar(pair(1,k)).eq.0) pair(1,k) = 'W'
	   if(iachar(pair(2,k)).eq.0) pair(2,k) = 'W'
	   if(iachar(pair(3,k)).eq.0) pair(3,k) = 'W'
C	write(*,*) infbp(1,k),infbp(2,k),iachar(pair(1,k)),iachar(pair(2,k))
c       if(INFBP(2,k).ne.0)then
         WRITE(11,431)INFBP(1,K),BASE(K),INFBP(2,K),bproll(K),bptilt(K),
     1   bptwst(K),bpslid(K),bpshft(K),bprise(K),pair(1,k),pair(2,k),
     2   pair(3,k),filprt
c       endif
 2015	CONTINUE
      CALL STATIS(bpTILT,LNDEX+1,AVTL,STDTL)
      CALL STATIS(bpROLL,LNDEX+1,AVRL,STDRL)
      CALL STATIS(bpTWST,LNDEX+1,AVTW,STDTW)
      CALL STATIS(bpshft,LNDEX+1,AVSX,STDSX)
      CALL STATIS(bpslid,LNDEX+1,AVSY,STDSY)
      CALL STATIS(bprise,LNDEX+1,AVSZ,STDSZ)
      WRITE(11,199)
      WRITE(11,110) avrl,avtl,avtw,avsy,avsx,avsz
      WRITE(11,111) stdrl,stdtl,stdtw,stdsy,stdsx,stdsz
      WRITE(11,199)
	if(npradd.gt.0) then
	  write(11,2017)
 2017	format(/'Basepair parameters for the multi-stranded regions'/)
	  do m=1,npradd
	    k=m+lndex+1
	    bpadd=badd1(m)//":"//badd2(m)
	    if(prtype(k).eq.'T') then
       WRITE(11,432)npadd1(m),bpadd,npadd2(m),bproll(K),bptilt(K),
     1      bptwst(K),bpslid(K),bpshft(K),bprise(K),pair(1,k),pair(2,k),
     2      pair(3,k),filprt
	 
	WRITE(11,434)infbp(1,npadd1(m)),base(npadd1(m)),infbp(2,
     1 npadd1(m)), bproll(npadd1(m)),bptilt(npadd1(m)),
     1      bptwst(npadd1(m)),bpslid(npadd1(m)),
     2 bpshft(npadd1(m)),bprise(npadd1(m)),pair(1,npadd1(m)),
     3 pair(2,npadd1(m)), pair(3,npadd1(m)),filprt
	 write(11,*) '  '
	    endif
	  enddo	  
          write(11,199)
	  do m=1,npradd
	    k=m+lndex+1
	    bpadd=badd1(m)//":"//badd2(m)
	    if(prtype(k).eq.'B') then
	      WRITE(11,433)npadd1(m),bpadd,npadd2(m),bproll(K),bptilt(K),
     1      bptwst(K),bpslid(K),bpshft(K),bprise(K),pair(1,k),pair(2,k),
     2      pair(3,k),filprt
	 
	WRITE(11,435)infbp(1,npadd1(m)),base(npadd1(m)),infbp(2,
     1 npadd1(m)), bproll(npadd1(m)),bptilt(npadd1(m)),
     1      bptwst(npadd1(m)),bpslid(npadd1(m)),
     2 bpshft(npadd1(m)),bprise(npadd1(m)),pair(1,npadd1(m)),
     3 pair(2,npadd1(m)), pair(3,npadd1(m)),filprt
c        WRITE(11,431)npadd1(m),bpadd,npadd2(m),bproll(m),bptilt(m),
c    1      bptwst(m),bpslid(m),bpshft(m),bprise(m),pair(1,m),pair(2,m),
c    2      pair(3,m),filprt
	      write(11,*) '  '
	    endif
	  enddo	  
          write(11,199)
	 endif
       END IF
        IF(NTIM.LE.1) THEN
          IF(ANSrnt .NE. 'Y')THEN
          WRITE(11,90)
        ELSE
 90   FORMAT(//'  Molecule has NOT been reoriented along GLOBAL axis')
       IF(ANORIN .EQ. 'C   ') THEN
          WRITE(11,92)
 92     FORMAT(//'  Molecule has been reoriented along best linear axis 
     1obtained from '/'  basepair centers.')
       ELSE IF(ANORIN.EQ.'O   ') THEN
          WRITE(11,91)
 91     FORMAT(//'  Molecule has been reoriented along best linear axis 
     1obtained from '/'  local helix ORIGINs.')
       ELSE
           WRITE(11,192)ANORIN
 192    FORMAT(//'  Molecule has been reoriented along best linear axis 
     1obtained from ',A4,' atoms.')
       END IF
        WRITE(11,93) DCRS1,DCRS2,DCRS3,ANGLE
        WRITE(11,94) DX1,DY1
       END IF
C
      END IF
900    CONTINUE
         WRITE(11,2)
  2   FORMAT(/2X,'GLOBAL Helical Parameters:',/20X,'Inclin.   Tip    Twi
     1st','     dx      dy      dz    ')
         WRITE(11,199)
 93   FORMAT(/2X,'Molecule rotated along(',3F8.4,') by',F8.3,' Deg.')
 94   FORMAT(2X,'and translated by',2F7.3,'A along X and Y directions.')
      DO 2003 K=1,LNDEX
C         IF(TWISTH(K) .LT. -0.001)TWISTG(K) = -TWISTG(K)
        if((anshlx.eq.'Y').or.(anshlx.eq.'y'))then
          do ih=1,infopr
            if((ihlx(1,ih).ne.0).and.(infbp(1,k).eq.ihlx(1,ih)))then
         WRITE(11,440)INFBP(1,K), BASE(K),INFBP(2,K),TILTG(K),ROLLG(K),
     1             TWISTG(K),DXG(K),DYG(K),HG(K),filprt
            endif
          enddo
        else
         WRITE(11,440)INFBP(1,K), BASE(K),INFBP(2,K),TILTG(K),ROLLG(K),
     1             TWISTG(K),DXG(K),DYG(K),HG(K),filprt
        endif

 2003 CONTINUE
      KK = LNDEX + 1
      CALL STATIS(TILTG,(LNDEX+1),AVTL,STDTL)
      CALL STATIS(ROLLG,(LNDEX+1),AVRL,STDRL)
      CALL STATIS(TWISTG,(LNDEX),AVTW,STDTW)
      CALL STATIS(DXG,(LNDEX+1),AVSX,STDSX)
      CALL STATIS(DYG,(LNDEX+1),AVSY,STDSY)
      CALL STATIS(HG,(LNDEX),AVSZ,STDSZ)
      WRITE(11,441) INFBP(1,KK),BASE(KK),INFBP(2,KK),TILTG(KK),ROLLG(KK)
     1,DXG(KK),DYG(KK),filprt
      WRITE(11,199)
      WRITE(11,110) AVTL,AVRL,AVTW,AVSX,AVSY,AVSZ
      WRITE(11,111) STDTL,STDRL,STDTW,STDSX,STDSY,STDSZ
      WRITE(11,199)
 950  CONTINUE
C-----------------------------------------------------------------------
C-- NOT PRINTED OUT FOR SINGLE STRAND CALCULATION.
 1000 CONTINUE
C-----------------------------------------------------------------------
c      IF (NTIM .EQ. 1.and.ansbpn.eq.'Y'.and.ansdb.eq.'Y')then
      IF (NTIM .EQ. 1.and.ansdb.eq.'Y')then
	WRITE(11,4)
c	write(*,*) 'Going to write axes and origins'
  4   FORMAT(//6X,' Local Helix Origins',15X,'Local Helix Axes'/)
        DO 2005 K=1,LNDEX
          ZS1 = ZAXISH(1,K) * 1.0
          ZS2 = ZAXISH(2,K) * 1.0
          ZS3 = ZAXISH(3,K) * 1.0
        IF(NTIM .EQ. 1)WRITE(11,200)(HORIG(K,M),M=1,3),ZS1,ZS2,ZS3
 2005	CONTINUE
  201  FORMAT(51X,3F8.3)
C
       IF(NTIM .EQ. 1)WRITE(11,96)
  96   FORMAT(//6X,' Base/B.P.  Centers',16X,'Base/B.P. Normals'/)
         DO 2006 K=1,LNDEX+1
           ZS1 = BPNORML(K,1) * 1.0
           ZS2 = BPNORML(K,2) * 1.0
           ZS3 = BPNORML(K,3) * 1.0
         IF(NTIM .EQ.1)WRITE(11,200)(BSCENTR(K,M),M=1,3),ZS1,ZS2,ZS3
 2006	CONTINUE
	endif
 200  FORMAT(3X,3F8.3,5X,3F10.4)
      INDM = LNDEX - 1
      DO 2007 K=1,LNDEX
         DO 2008 KK=1,LNDEX
           SUM = 0.0
           DO 2009 KJ=1,3
              SUM = SUM + ZAXISH(KJ,K) * ZAXISH(KJ,KK)
 2009	    CONTINUE
           ANG = ACOS(SUM) * CONV
           ANGST(K,KK) = ANG
 2008	  CONTINUE
 2007	CONTINUE
c        WRITE(11,97)
c        DO 2010 KKK=1,LNDEX
c           WRITE(11,11) KKK,(ANGST(KKK,K),K=1,LNDEX)
c 2010   CONTINUE
	WRITE(11,99) ANGST(LNDEX,1)
  99	FORMAT(/3X,'END-TO-END BENDING ANGLE:',F7.2/)
  97    FORMAT(//3X,'Angle between successive Local Helix Axes:'/)
  11    FORMAT(2X,I2,'.',10F7.2,9(/5X,10F7.2))
C
         DO 2011 K=1,LNDEX + 1
           DO 2012 KK=1,LNDEX + 1
             SUM = 0.0
             DO 2013 KJ = 1,3
               SUM = SUM + BPNORML(K,KJ)*BPNORML(KK,KJ)
 2013        CONTINUE
             ANGST(K,KK) = ACOSF(SUM)
 2012	    CONTINUE
 2011	 CONTINUE
  98    FORMAT(/3X,'Angle between successive Base/Base-pair Normals:'/)
c        WRITE(11,98)
        DO 2014 KKK=1,LNDEX + 1
c           WRITE(11,11) KKK,(ANGST(KKK,K),K=1,LNDEX+1)
 2014	CONTINUE
      RETURN
      END
C----------------------------------------------------------------------
      SUBROUTINE SLIDCAL
	parameter ( nrs = 9999 )
	include 'parameters.h'
      DOUBLE PRECISION XAXIS,YAXIS,ZAXIS,XM,YM,ZM1,ZM2,ZM3,Y1,Y2,Y3,
     1   Y12,Y22,Y32,X11,X12,X13,BMN,TL2,RL2,WTILT,WROLL,Y11,Y13
      DIMENSION N(3),M(3),A(3,3)
C 
      SLIDEX = SLX(LNDEX)
      SLIDEY = SLY(LNDEX)
      SLIDEZ = SLZ(LNDEX)
      TILT = DIETLT(LNDEX)
      ROLL = DIEROL(LNDEX)
      TWIST = TWISTL(LNDEX)
      CONV = 3.141592654/180.000
      ROLL = ROLL * CONV
      TILT = TILT * CONV
      TWIST = TWIST * CONV
      CY = COS(ROLL)
      CX = COS(TILT)
      SY = SIN(-ROLL)
      SX = SIN(TILT)
      CT = COS(TWIST)
      ST = SIN(TWIST)
      A(1,1) = 2.0 * CX * ST
      A(1,2) = 0.0
      A(1,3) = 2.0 * SX
      A(2,1) = 0.0
      A(2,2) = -2.0 * CY * ST
      A(2,3) = 2.0 * SY
      A(3,1) = -2.0 * CY * SX * ST
      A(3,2) = 2.0 * CX * SY * ST
      A(3,3) = CX * CY * (1.0 + CT)
      CALL MINV(A,3,DXX,N,M)
      DX = 1.0 + CT + SX * SX * ( 1.0 - CT)
      DX = SQRT(2.0 * DX)
      DY = 1.0 + CT + SY * SY * (1.0 - CT)
      DY = SQRT(2.0 * DY)
      B1 = SLIDEY * DY
      B2 = SLIDEX * DX
      B3 = (SLIDEZ * DX * DY)/2
      D1 = A(1,1) * B1 + A(1,2) * B2 + A(1,3) * B3
      D2 = A(2,1) * B1 + A(2,2) * B2 + A(2,3) * B3
      D3 = A(3,1) * B1 + A(3,2) * B2 + A(3,3) * B3
      DXLOC(LNDEX) = D1
      DYLOC(LNDEX) = D2
      DZLOC(LNDEX) = D3
      RETURN
      END
C
C***********************************************************************
      SUBROUTINE MINV (A,N,D,L,M)
C***********************************************************************
C     ----- STANDARD IBM MATRIX INVERSION ROUTINE -----
C
      DIMENSION A(1),L(1),M(1)
C
C     ----- SEARCH FOR LARGEST ELEMENT -----
C
      D = 1.E0
      NK = -N
      DO 80 K = 1,N
      NK = NK+N
      L(K) = K
      M(K) = K
      KK = NK+K
      BIGA = A(KK)
      DO 20 J = K,N
      IZ = N*(J-1)
      DO 20 I = K,N
      IJ = IZ+I
c   10 IF( ABS(BIGA)- ABS(A(IJ))) 15,20,20
   10 continue
	IF( ( ABS(BIGA)- ABS(A(IJ)) ) .lt.0.0) then
   15 BIGA = A(IJ)
      L(K) = I
      M(K) = J
      endif
   20 CONTINUE
C
C     ----- INTERCHANGE ROWS -----
C
      J = L(K)
c      IF(J-K) 35,35,25
      IF( (J-K).gt.0) then
   25 KI = K-N
      DO 30 I = 1,N
      KI = KI+N
      HOLD = -A(KI)
      JI = KI-K+J
      A(KI) = A(JI)
   30 A(JI) = HOLD
C
C     ----- INTERCHANGE COLUMNS -----
C
      endif
   35 I = M(K)
c      IF(I-K) 45,45,38
      IF( (I-K) .gt.0) then
   38 JP = N*(I-1)
      DO 40 J = 1,N
      JK = NK+J
      JI = JP+J
      HOLD = -A(JK)
      A(JK) = A(JI)
   40 A(JI) = HOLD
C
C     ----- DIVIDE COLUMN BY MINUS PIVOT -----
C
c   45 IF(BIGA) 48,46,48
c   46 D = 0.E0
       endif
   45 IF(BIGA .eq.0) then
   46 D = 0.E0
      GO TO 150
      endif
   48 DO 55 I = 1,N
c      IF(I-K) 50,55,50
      IF( (I-K) .ne. 0) then         ! 50,55,50
   50 IK = NK+I
      A(IK) = A(IK)/(-BIGA)
      endif
   55 CONTINUE
C
C     ----- REDUCE MATRIX -----
C
      DO 65 I = 1,N
      IK = NK+I
      HOLD = A(IK)
      IJ = I-N
      DO 65 J = 1,N
      IJ = IJ+N
c      IF(I-K) 60,65,60
      IF( (I-K) .ne.0 ) then
c   60 IF(J-K) 62,65,62
   60 IF( (J-K) .ne. 0) then
   62 KJ = IJ-I+K
      A(IJ) = HOLD*A(KJ)+A(IJ)
      endif
      endif
   65 CONTINUE
C
C     ----- DIVIDE ROW BY PIVOT -----
C
      KJ = K-N
      DO 75 J = 1,N
      KJ = KJ+N
c      IF(J-K) 70,75,70
      IF( (J-K) .ne.0) then
   70 A(KJ) = A(KJ)/BIGA
      endif
   75 CONTINUE
C
C     ----- PRODUCT OF PIVOTS -----
C
      D = D*BIGA
C
C     ----- REPLACE PIVOT BY RECIPROCAL -----
C
      A(KK) = 1.E0/BIGA
   80 CONTINUE
C
C     ----- FINAL ROW AND COLUMN INTERCHANGE -----
C
      K = N
  100 K = (K-1)
c      IF(K) 150,150,105
      IF( K.gt.0) then   ! 150,150,105
  105 I = L(K)
c      IF(I-K) 120,120,108
      IF( (I-K) .gt.0 ) then    ! 120,120,108
  108 JQ = N*(K-1)
      JR = N*(I-1)
      DO 110 J = 1,N
      JK = JQ+J
      HOLD = A(JK)
      JI = JR+J
      A(JK) = -A(JI)
  110 A(JI) = HOLD
      endif
  120 J = M(K)
c      IF(J-K) 100,100,125
      IF( (J-K) .le.0 ) goto 100   ! 100,100,125
  125 KI = K-N
      DO 130 I = 1,N
      KI = KI+N
      HOLD = A(KI)
      JI = KI-K+J
      A(KI) = -A(JI)
  130 A(JI) = HOLD
      GOTO 100
      endif
  150 RETURN
      END
C
      SUBROUTINE SINGAXS(IP,N,NATL,NC8,X,L,XAXIS,Y,EL,EM,EN,IBASE,SDVX)
	parameter ( natm = 999999 )
	parameter ( nrs = 9999 )
	include 'coordinates.h'
      DOUBLE PRECISION XAX,XAXIS,EL,EM,EN,Y,TEMPX,TEMPY
      DIMENSION IP(4),IN(4),XAX(3),YAXIS(2,3),ORIG(3),X(3,15),Y(2,3),
     1          XAXIS(2,3),N(30),DVA(30),TEMPX(3),TEMPY(3),IBASE(10)
     2         ,R(3,3)
C
      L2 = L * 2
      NADD = 0
      DO 2000 I=1,NATL
          IF(IBASE(I).EQ.IP(L2)) NC1P = I
 2000 CONTINUE
      DO 2001 KJ=1,3
         XAX(KJ) = X(KJ,NC8) - X(KJ,NC1P)
 2001 CONTINUE
      CALL NORMAL(XAX(1),XAX(2),XAX(3),XAXIS(L,1),XAXIS(L,2),XAXIS(L,3))
	if(natl.ge.5) then
      CALL PLANEB(X,NATL,N,1,DVA,EL,EM,EN,DETR,SDVX)
	else
	  write(6,801) natl
	  stop
	endif
      CALL CROSS(EL,EM,EN,XAXIS(L,1),XAXIS(L,2),XAXIS(L,3),Y(L,1),Y(L,2)
     1          ,Y(L,3))
      RETURN
 801	format('BAD Selection of BASE, only',I3,' atoms found')
      END
C
      SUBROUTINE ROTVEC(RMAT,A1,A2,A3,B)
      DOUBLE PRECISION A1,A2,A3,B,DRMAT(3,3)
      DIMENSION RMAT(3,3),B(3)
 
C       CONVERT MATRIX TO DOUBLE PRECISION
 
        DO 2002 I=1,3
          DO 2003 J=1,3
             DRMAT(I,J)=DBLE(RMAT(I,J))
 2003	   CONTINUE
 2002	 CONTINUE   
        
      DO 2004 I=1,3
         B(I) = DRMAT(I,1) * A1 + DRMAT(I,2) * A2 + DRMAT(I,3) * A3
 2004 CONTINUE
 
      RETURN
      END
C
      SUBROUTINE STATIS(VAR,NVAR,AVG,STD)
      DIMENSION VAR(1)
C
      IF(NVAR.GE.2) THEN
         SUM = 0.0
         NPOINT = 0
         DO 10 IV = 1,NVAR
           IF(VAR(IV).NE.0.000000) THEN
             SUM = SUM + VAR(IV)
   	      NPOINT = NPOINT + 1
   	    END IF
10       CONTINUE
         AVG = SUM/NPOINT
         SUM = 0.0
         DO 20 IV = 1, NVAR
	    IF(VAR(IV).NE.0.000000) THEN
             SUM = SUM + VAR(IV)**2
  	    END IF
20      CONTINUE
         SUM = SUM/NPOINT
         STD = SUM-AVG**2
         IF(STD .LT. 0.0) STD = 0.0
     	  IF(NPOINT.GE.1) STD = SQRTF(STD)
      ELSE
           AVG = VAR(1)
           STD = 0.0
      END IF
      RETURN
      END
C
      SUBROUTINE ROTMAT(EL,EM,EN,THETA,R)
      DIMENSION R(3,3)
C
      CON=3.14195/180.0
      THETR=THETA*CON
      CST= COS(THETR)
      SNT= SIN(THETR)
      CSTF=1.0-CST
      R(1,1)=EL*EL*CSTF+CST
      R(2,2)=EM*EM*CSTF+CST
      R(3,3)=EN*EN*CSTF+CST
      AA=EL*EM*CSTF
      BB=EN*SNT
      R(1,2)=AA-BB
      R(2,1)=AA+BB
      AA=EL*EN*CSTF
      BB=EM*SNT
      R(1,3)=AA+BB
      R(3,1)=AA-BB
      AA=EM*EN*CSTF
      BB=EL*SNT
      R(2,3)=AA-BB
      R(3,2)=AA+BB
      RETURN
      END
C
        SUBROUTINE MATMUL(R,A1,A2,A3,B)
        DIMENSION B(3),R(3,3)
C
        DO 2000 I=1,3
            B(I) = R(I,1) * A1 + R(I,2) * A2 + R(I,3) * A3
 2000	 CONTINUE
        RETURN
        END
C
        SUBROUTINE POINT(NPTS,DX,DY,CENTRS)
	parameter ( nrs = 9999 )
        DIMENSION CENTRS(nrs,3)
C
        DX = 0.0
        DY = 0.0
        DO 2000 KJ=1,NPTS
           DX = DX + CENTRS(KJ,1)
           DY = DY + CENTRS(KJ,2)
 2000	 CONTINUE
        DX = DX/NPTS
        DY = DY/NPTS
        RETURN
        END
C
      SUBROUTINE REWRIT(FILENM,TITLE,KKTI)
	parameter ( natm = 999999 )
	parameter ( nrs = 9999 )
	include 'coordinates.h'
      CHARACTER*80 TITLE(20)
      CHARACTER*40 FILENM,FILNEW
     
C
        KKKK=INDEX(FILENM(1:40),':')
        IF (KKKK .EQ. 0) THEN
         KKKK=1
        ELSE
         KKKK=KKKK+1
        ENDIF
C
        KKKKK=INDEX(FILENM(KKKK:40),'.')
        KKKKK=KKKKK+KKKK
C
        FILNEW=FILENM(KKKK:KKKKK-1)//'coor'
C
      OPEN(UNIT=7,FILE=FILNEW)
C
      DO 2000 I=1,KKTI
               WRITE(7,'(A40)') TITLE(I)
 2000 CONTINUE
      K1 = 1
      NEND = IENS(K1)
      DO 2001 K=1,NAT
             IF(K.GT.NEND) THEN
                   K1 = K1 + 1
                   NEND = IENS(K1)
             END IF
             WRITE(7,101) ATNM(K),SECNM(K1),K1,(CR(K,M),M=1,3)
 2001 CONTINUE
  101 FORMAT('ATOM',8X,A4,1X,A3,2X,I4,4X,3F8.3,18x,'MODL')
C       WRITE(6,10)FILNEW
       WRITE(8,10)FILNEW
  10  FORMAT(1X,'Reoriented atomic coordinates written into:'/1X,A80)
       CLOSE(7)
       RETURN
       END
C
C-----------------------------------------------------------------------
       SUBROUTINE TORSION(ANSDB,ANSPDB,filprt)
	parameter ( natm = 999999 )
	parameter ( nrs = 9999 )
	include 'coordinates.h'
	dimension ac41(nrs,3),ap2(nrs,3)
       DIMENSION COOR1(3),COOR2(3),N(5),COOR3(3),NN(natm,6),ALPHA(nrs),
     1        XYZ1(3),XYZ2(3),XYZ3(3),XYZ4(3),BETA(nrs),GAMA(nrs),
     2        DELTA(nrs),EPSILN(nrs),ZETA(nrs),CHI(nrs),ANEU0(nrs),
     3        ANEU1(nrs),ANEU2(nrs),ANEU3(nrs),ANEU4(nrs),ANEUMAX(nrs),
     4        PHASE(nrs),XYZ(natm,3),EE(21),ISTSN(nrs),IENSN(nrs),
     5        pseta(nrs),pstheta(nrs),SECNMN(nrs)
	CHARACTER*1 ANSDB,ANSPDB,ANSAMB
       CHARACTER*12 EE,CLASS,filprt*20
       CHARACTER*3 SECNMN
       CHARACTER*4 AN1,AN2,AN3,AN4,AC5,AO5,AC4,AC3,AO3,AP,AO4,
     1              AC2,AC1,AO1
       DOUBLE PRECISION DIST,DIS,EL
       EQUIVALENCE (XYZ(1,1),CR(1,1))
       logical pfound
C
       DATA EE/'2*EXO-3*ENDO','3*ENDO','3*ENDO-4*EXO','4*EXO',
     1   '4*EXO-O*ENDO','O*ENDO','O*ENDO-1*EXO','1*EXO','1*EXO-2*ENDO',
     2   '2*ENDO', '2*ENDO-3*EXO','3*EXO','3*EXO-4*ENDO','4*ENDO',
     3   '4*ENDO-O*EXO','O*EXO','O*EXO-1*ENDO','1*ENDO','1*ENDO-2*EXO',
     4   '2*EXO','2*EXO-3*ENDO'/
 
 1111   FORMAT(A20)
C
C         FINDING THE CONNECTIVITY MATRIX OF THE ATOMS
C
      DATA AP,AO5,AC5,AC4,AC3,AO3/'P   ','O5'' ','C5'' ','C4'' ','C3'' '
     1                           ,'O3'' '/
        DATA AC2,AC1,AO1,AO4/'C2'' ','C1'' ','O1'' ','O4'' '/
        CONV = 3.14159/180.0
        DO 2000 KS=1,NSEG
             ALPHA(KS) = 0.0
             BETA(KS) = 0.0
             GAMA(KS) = 0.0
             DELTA(KS) = 0.0
             EPSILN(KS) = 0.0
             ZETA(KS) = 0.0
             CHI(KS) = 0.0      
             pseta(KS) = 0.0      
             pstheta(KS) = 0.0      
             ANEU0(KS) = 0.0
             ANEU1(KS) = 0.0
             ANEU2(KS) = 0.0
             ANEU3(KS) = 0.0
             ANEU4(KS) = 0.0
 2000	 CONTINUE
	IF(ANSPDB.NE.'Y') THEN
C	   WRITE(6,1001)
 1001	FORMAT(/'Is the data file in AMBER format? [Y]',$)
C	   READ(5,'(A1)') ANSAMB
C	   IF(ANSAMB.NE.'n'.AND.ANSAMB.NE.'N') ANSAMB = 'Y'
	   ansamb = 'N'
	   IF(ANSAMB.EQ.'Y') THEN
	      KSEG = 0
	      DO 2101 KRES=1,NSEG
	         IF(SECNM(KRES).EQ.'HB '.OR.SECNM(KRES).EQ.'POM') THEN
	            KSEG = KSEG + 1
	            IENSN(KSEG) = IENS(KRES+1)
	            ISTSN(KSEG) = ISTS(KRES)
	            SECNMN(KSEG) = SECNM(KRES+1)
	         END IF
	         IF(SECNM(KRES).EQ.'HE ') THEN
	            SECNMN(KSEG) = SECNM(KRES-1)
	            IENSN(KSEG) = IENS(KRES)
	         END IF
 2101	      CONTINUE
	      DO 2102 KRES=1,NSEG
	        ISTS(KRES) = 0
	        IENS(KRES) = 0
	        SECNM(KRES) = '   '
 2102	      CONTINUE
	      DO 2103 KRES=1,KSEG
	        ISTS(KRES) = ISTSN(KRES)
	        IENS(KRES) = IENSN(KRES)
	        SECNM(KRES) = SECNMN(KRES)
 2103	     CONTINUE
	     NSEG = KSEG
	   END IF
	 END IF
C
           DO 2001 I=1,NAT
              DO 2002 J=1,6
                 NN(I,J) = 0
 2002	       CONTINUE
 2001 	    CONTINUE
C
        DO 10 I=1,NAT
            KJ=0
            DO 2003 J=1,(I-1)
              DO 2004 KN=1,6
                IF(I.EQ.NN(J,KN)) THEN
                   KJ = KJ + 1
                   NN(I,KJ) = J
                END IF
 2004	       CONTINUE
 2003	     CONTINUE
            DO 2005 J=(I+1),(I+70)
	         IF(J.GT.NAT) GO TO 10
                   DIST=0.0
                   DO 2006 K=1,3
                       DIST=DIST+(XYZ(I,K)-XYZ(J,K))*(XYZ(I,K)-XYZ(J,K))
 2006	            CONTINUE
                   IF(DIST .LT. 3.5  .AND. KJ .LT. 6) THEN
                        KJ=KJ+1
                        NN(I,KJ)=J
                   END IF
 2005	      CONTINUE
  10    CONTINUE
        DO 100 I=1,NAT
            AN1 = ATNM(I)
           DO 200 K=1,6
              N1=NN(I,K)
              IF(N1 .EQ. 0)GO TO 200
              DO 2007 KS=1,NSEG
                 IF(N1.GE.ISTS(KS).AND.N1.LT.IENS(KS)) THEN
                    IUNIT = KS
                 END IF
 2007	       CONTINUE
            AN2 = ATNM(N1)
              DO 300 L=1,6
                  N2=NN(N1,l)
                 IF(N2.LE.N1.OR.N2.EQ.I)GO TO 300
                    AN3 = ATNM(N2)
                  DO 400 M=1,6
                       N3=NN(N2,M)
                       IF(N3.EQ.0.OR.N3.EQ.N1)GO TO 400
                          AN4 = ATNM(N3)
                       DO 2008 IJ=1,3
                           XYZ1(IJ)=XYZ(I,IJ)
                           XYZ2(IJ)=XYZ(N1,IJ)
                           XYZ3(IJ)=XYZ(N2,IJ)
                           XYZ4(IJ)=XYZ(N3,IJ)
 2008	                CONTINUE
      IF(AN1.EQ.AO5.AND.AN2.EQ.AC5.AND.AN3.EQ.AC4.AND.AN4.EQ.AC3) THEN
              CALL TORAN(XYZ1,XYZ2,XYZ3,XYZ4,TAU)
              GAMA(IUNIT) = TAU
      END IF
      IF(AN1.EQ.AC3.AND.AN2.EQ.AC4.AND.AN3.EQ.AC5.AND.AN4.EQ.AO5) THEN
              CALL TORAN(XYZ1,XYZ2,XYZ3,XYZ4,TAU)
              GAMA(IUNIT) = TAU
      END IF
      IF(AN1.EQ.AP.AND.AN2.EQ.AO5.AND.AN3.EQ.AC5.AND.AN4.EQ.AC4) THEN
              CALL TORAN(XYZ1,XYZ2,XYZ3,XYZ4,TAU)
              BETA(IUNIT) = TAU
      END IF
      IF(AN1.EQ.AC4.AND.AN2.EQ.AC5.AND.AN3.EQ.AO5.AND.AN4.EQ.AP) THEN
              CALL TORAN(XYZ1,XYZ2,XYZ3,XYZ4,TAU)
              BETA(IUNIT) = TAU
      END IF
      IF(AN1.EQ.AC5.AND.AN2.EQ.AC4.AND.AN3.EQ.AC3.AND.AN4.EQ.AO3) THEN
              CALL TORAN(XYZ1,XYZ2,XYZ3,XYZ4,TAU)
              DELTA(IUNIT) = TAU
      END IF
      IF(AN1.EQ.AO3.AND.AN2.EQ.AC3.AND.AN3.EQ.AC4.AND.AN4.EQ.AC5) THEN
              CALL TORAN(XYZ1,XYZ2,XYZ3,XYZ4,TAU)
              DELTA(IUNIT) = TAU
      END IF
      IF(AN1.EQ.AC4.AND.AN2.EQ.AC3.AND.AN3.EQ.AO3.AND.AN4.EQ.AP) THEN
              CALL TORAN(XYZ1,XYZ2,XYZ3,XYZ4,TAU)
              EPSILN(IUNIT) = TAU
      END IF
      IF(AN1.EQ.AP.AND.AN2.EQ.AO3.AND.AN3.EQ.AC3.AND.AN4.EQ.AC4) THEN
              CALL TORAN(XYZ1,XYZ2,XYZ3,XYZ4,TAU)
              EPSILN(IUNIT) = TAU
      END IF
      IF(AN1.EQ.AC3.AND.AN2.EQ.AO3.AND.AN3.EQ.AP.AND.AN4.EQ.AO5) THEN
              CALL TORAN(XYZ1,XYZ2,XYZ3,XYZ4,TAU)
              ZETA(IUNIT) = TAU
      END IF
      IF(AN1.EQ.AO5.AND.AN2.EQ.AP.AND.AN3.EQ.AO3.AND.AN4.EQ.AC3) THEN
              CALL TORAN(XYZ1,XYZ2,XYZ3,XYZ4,TAU)
              ZETA(IUNIT) = TAU
      END IF
      IF(AN1.EQ.AO3.AND.AN2.EQ.AP.AND.AN3.EQ.AO5.AND.AN4.EQ.AC5) THEN
              CALL TORAN(XYZ1,XYZ2,XYZ3,XYZ4,TAU)
              ALPHA(IUNIT) = TAU
      END IF
      IF(AN1.EQ.AC5.AND.AN2.EQ.AO5.AND.AN3.EQ.AP.AND.AN4.EQ.AO3) THEN
              CALL TORAN(XYZ1,XYZ2,XYZ3,XYZ4,TAU)
              ALPHA(IUNIT) = TAU
      END IF
      IF((AN1.EQ.'O4'' '.OR.AN1.EQ.'O1'' ').AND.AN2.EQ.'C1'' ') THEN
             IF(AN3.EQ.'N9  '.AND.AN4.EQ.'C4  ') THEN
                 CALL TORAN(XYZ1,XYZ2,XYZ3,XYZ4,TAU)
                 CHI(IUNIT) = TAU
             END IF
             IF(AN3.EQ.'N1  '.AND.AN4.EQ.'C2  ') THEN
                 CALL TORAN(XYZ1,XYZ2,XYZ3,XYZ4,TAU)
                 CHI(IUNIT) = TAU
             END IF
      ELSE IF((AN4.EQ.'O4'' '.OR.AN4.EQ.'O1'' ').AND.AN3.EQ.'C1'' ')THEN
             IF(AN2.EQ.'N9  '.AND.AN1.EQ.'C4  ') THEN
                 CALL TORAN(XYZ1,XYZ2,XYZ3,XYZ4,TAU)
                 CHI(IUNIT) = TAU
             END IF
             IF(AN2.EQ.'N1  '.AND.AN1.EQ.'C2  ') THEN
                 CALL TORAN(XYZ1,XYZ2,XYZ3,XYZ4,TAU)
                 CHI(IUNIT) = TAU
             END IF
      END IF
        IF(AN1.EQ.AC4.AND.AN2.EQ.AO1.AND.AN3.EQ.AC1.AND.AN4.EQ.AC2)THEN
              CALL TORAN(XYZ1,XYZ2,XYZ3,XYZ4,TAU)
              ANEU0(IUNIT) = TAU
        END IF
        IF(AN1.EQ.AC2.AND.AN2.EQ.AC1.AND.AN3.EQ.AO1.AND.AN4.EQ.AC4) THEN
              CALL TORAN(XYZ1,XYZ2,XYZ3,XYZ4,TAU)
              ANEU0(IUNIT) = TAU
        END IF
        IF(AN1.EQ.AC4.AND.AN2.EQ.AO4.AND.AN3.EQ.AC1.AND.AN4.EQ.AC2)THEN
              CALL TORAN(XYZ1,XYZ2,XYZ3,XYZ4,TAU)
              ANEU0(IUNIT) = TAU
        END IF
        IF(AN1.EQ.AC2.AND.AN2.EQ.AC1.AND.AN3.EQ.AO4.AND.AN4.EQ.AC4) THEN
              CALL TORAN(XYZ1,XYZ2,XYZ3,XYZ4,TAU)
              ANEU0(IUNIT) = TAU
        END IF
        IF(AN1.EQ.AC3.AND.AN2.EQ.AC2.AND.AN3.EQ.AC1.AND.AN4.EQ.AO1)THEN
              CALL TORAN(XYZ1,XYZ2,XYZ3,XYZ4,TAU)
              ANEU1(IUNIT) = TAU
        END IF
        IF(AN1.EQ.AO1.AND.AN2.EQ.AC1.AND.AN3.EQ.AC2.AND.AN4.EQ.AC3) THEN
              CALL TORAN(XYZ1,XYZ2,XYZ3,XYZ4,TAU)
              ANEU1(IUNIT) = TAU
        END IF
        IF(AN1.EQ.AC3.AND.AN2.EQ.AC2.AND.AN3.EQ.AC1.AND.AN4.EQ.AO4)THEN
              CALL TORAN(XYZ1,XYZ2,XYZ3,XYZ4,TAU)
              ANEU1(IUNIT) = TAU
        END IF
        IF(AN1.EQ.AO1.AND.AN2.EQ.AC1.AND.AN3.EQ.AC2.AND.AN4.EQ.AC3)THEN
           CALL TORAN(XYZ1,XYZ2,XYZ3,XYZ4,TAU)
           ANEU1(IUNIT) = TAU
        END IF
        IF(AN1.EQ.AO4.AND.AN2.EQ.AC1.AND.AN3.EQ.AC2.AND.AN4.EQ.AC3)THEN
          CALL TORAN(XYZ1,XYZ2,XYZ3,XYZ4,TAU)
          ANEU1(IUNIT) = TAU
        END IF
        IF(AN1.EQ.AC4.AND.AN2.EQ.AC3.AND.AN3.EQ.AC2.AND.AN4.EQ.AC1)THEN
           CALL TORAN(XYZ1,XYZ2,XYZ3,XYZ4,TAU)
           ANEU2(IUNIT) = TAU
        END IF
        IF(AN1.EQ.AC1.AND.AN2.EQ.AC2.AND.AN3.EQ.AC3.AND.AN4.EQ.AC4) THEN
              CALL TORAN(XYZ1,XYZ2,XYZ3,XYZ4,TAU)
              ANEU2(IUNIT) = TAU
        END IF
        IF(AN1.EQ.AO4.AND.AN2.EQ.AC4.AND.AN3.EQ.AC3.AND.AN4.EQ.AC2)THEN
           CALL TORAN(XYZ1,XYZ2,XYZ3,XYZ4,TAU)
           ANEU3(IUNIT) = TAU
        END IF
        IF(AN1.EQ.AC2.AND.AN2.EQ.AC3.AND.AN3.EQ.AC4.AND.AN4.EQ.AO4) THEN
              CALL TORAN(XYZ1,XYZ2,XYZ3,XYZ4,TAU)
              ANEU3(IUNIT) = TAU
        END IF
        IF(AN1.EQ.AO1.AND.AN2.EQ.AC4.AND.AN3.EQ.AC3.AND.AN4.EQ.AC2)THEN
           CALL TORAN(XYZ1,XYZ2,XYZ3,XYZ4,TAU)
           ANEU3(IUNIT) = TAU
        END IF
        IF(AN1.EQ.AC2.AND.AN2.EQ.AC3.AND.AN3.EQ.AC4.AND.AN4.EQ.AO1) THEN
              CALL TORAN(XYZ1,XYZ2,XYZ3,XYZ4,TAU)
              ANEU3(IUNIT) = TAU
        END IF
        IF(AN1.EQ.AC3.AND.AN2.EQ.AC4.AND.AN3.EQ.AO4.AND.AN4.EQ.AC1)THEN
           CALL TORAN(XYZ1,XYZ2,XYZ3,XYZ4,TAU)
           ANEU4(IUNIT) = TAU
        END IF
        IF(AN1.EQ.AC1.AND.AN2.EQ.AO4.AND.AN3.EQ.AC4.AND.AN4.EQ.AC3) THEN
              CALL TORAN(XYZ1,XYZ2,XYZ3,XYZ4,TAU)
              ANEU4(IUNIT) = TAU
        END IF
        IF(AN1.EQ.AC3.AND.AN2.EQ.AC4.AND.AN3.EQ.AO1.AND.AN4.EQ.AC1)THEN
           CALL TORAN(XYZ1,XYZ2,XYZ3,XYZ4,TAU)
           ANEU4(IUNIT) = TAU
        END IF
        IF(AN1.EQ.AC1.AND.AN2.EQ.AO1.AND.AN3.EQ.AC4.AND.AN4.EQ.AC3) THEN
              CALL TORAN(XYZ1,XYZ2,XYZ3,XYZ4,TAU)
              ANEU4(IUNIT) = TAU
        END IF
C
 400    CONTINUE
 300    CONTINUE
 200    CONTINUE
 100    CONTINUE
C
        nc4ser=0
	npser=0
	nresser=1
	IUNIT=2
	do nresser=1,nseg
	PFOUND=.FALSE.
          do i=ists(nresser),iens(nresser)
	  if(atnm(i).eq.'C4'' ') then
	    do l=1,3
	      ac41(nresser,l)=cr(i,l)
	    enddo
	   elseif(atnm(i).eq.'P   ') then
	     pfound=.TRUE.
	     do l=1,3
	       ap2(nresser,l)=cr(i,l)
	     enddo
	    endif
	  enddo
c	  write(*,*) nresser,ac41(nresser,1),ac41(nresser,2),ac41(nresser,3)
c	  write(*,*)'P',nresser,ap2(nresser,1),ap2(nresser,2),ap2(nresser,3)
        enddo
	do i=2,nseg-1
          do IJ=1,3
            XYZ1(IJ)=ac41(i-1,IJ)
            XYZ2(IJ)=ap2((i),IJ)
            XYZ3(IJ)=ac41((i),IJ)
            XYZ4(IJ)=ap2((i+1),IJ)
          enddo 
          dd1=distn(xyz1,xyz2)
          dd2=distn(xyz2,xyz3)
          dd3=distn(xyz3,xyz4)
c		write(*,*) iunit, ' ', dd1, ' ', dd2, ' ',dd3
          if(dd1.le.4.5.and.dd2.le.4.5.and.dd3.le.4.5) then
             CALL TORAN(XYZ1,XYZ2,XYZ3,XYZ4,TAU)
c	      write(*,*) 'Distances',dd1,dd2,dd3, 'Torsion=',iunit,tau ! not w
             pseta(IUNIT) = TAU
          else 
             pseta(iunit)=0.0
c	WRITE(*,*) "eta",iunit,TAU  ! not w
         endif
         IUNIT=IUNIT+1
         enddo
         IUNIT=2
         do kl=2,nseg-1
	do IJ=1,3
           XYZ1(IJ)=ap2(kl,IJ)
           XYZ2(IJ)=ac41(kl,IJ)
           XYZ3(IJ)=ap2((kl+1),IJ)
           XYZ4(IJ)=ac41((kl+1),IJ)
	enddo
        dd1=distn(xyz1,xyz2)
	dd2=distn(xyz2,xyz3)
	dd3=distn(xyz3,xyz4)
	if(dd1.le.4.5.and.dd2.le.4.5.and.dd3.le.4.5) then
	  CALL TORAN(XYZ1,XYZ2,XYZ3,XYZ4,TAU)
          pstheta(IUNIT) = TAU
c			WRITE(*,*) "theta",iunit,TAU
          else
	    pstheta(iunit)=0.0
	endif
	 IUNIT=IUNIT+1
	 enddo
                
        WRITE(11,110)
 110    FORMAT(//,2X,'BACKBONE AND GLYCOSIDIC TORSION ANGLES'/)
        WRITE(11,111)
 111    FORMAT(12X,' P-O5''  O5''-C5''  C5''-C4''  C4''-C3''  C3''-O3''
     1 O3''- P    C1''-N    Pseudo Torsions')
        WRITE(11,222)
 222    FORMAT(13X,'ALPHA     BETA    GAMMA    DELTA     EPS      ZETA
     1  CHI      ETA     THETA'/'  Eta & Theta are calculated following'
     2  ,' C.M. Duarte & A.M. Pyle (1998) J. Mol. Biol. 284,1465-1478')
        DO 2009 KS=1,NSEG
          WRITE(11,112) ks,secnm(ks)(1:1),ALPHA(KS),BETA(KS),GAMA(KS),
     1      DELTA(KS), EPSILN(KS),ZETA(KS),CHI(KS),pseta(ks),
     2      pstheta(ks),filprt
112       format('TR',i5,1x,a1,9F9.1,1x,a20) 
C
         IF(ALPHA(KS) .LT. 0.00)ALPHA(KS) = ALPHA(KS) + 360.0
         IF(BETA(KS) .LT. 0.00)BETA(KS) = BETA(KS) + 360.0
         IF(GAMA(KS) .LT. 0.00)GAMA(KS) = GAMA(KS) + 360.0
         IF(DELTA(KS) .LT. 0.00)DELTA(KS) = DELTA(KS) + 360.0
         IF(EPSILN(KS) .LT. 0.00)EPSILN(KS) = EPSILN(KS) + 360.0
         IF(ZETA(KS) .LT. 0.00)ZETA(KS) = ZETA(KS) + 360.0
         IF(CHI(KS) .LT. 0.00)CHI(KS) = CHI(KS) + 360.0
         IF(pseta(KS) .LT. 0.00)pseta(KS) = pseta(KS) + 360.0
         IF(pstheta(KS) .LT. 0.00)pstheta(KS) = pstheta(KS) + 360.0
 2009	CONTINUE
        write(11,199)
199     format(' The Pseudotorsions, namely ETA and THETA, are from',
     1 ' definition of C.M.Duarte and A.M.Pyle (1998) J.Mol.Biol.284,',
     2 '1465.')
       NSEG1 = NSEG - 1
      IF(ANSDB .EQ. 'Y') NSEG1 = NSEG - 2
         CALL STATIS(ALPHA,NSEG,AVALF,STDALF)
         CALL STATIS(BETA,NSEG,AVBET,STDBET)
         CALL STATIS(GAMA,NSEG,AVGAM,STDGAM)
         CALL STATIS(DELTA,NSEG,AVDEL,STDDEL)
         CALL STATIS(EPSILN,NSEG,AVEPS,STDEPS)
         CALL STATIS(ZETA,NSEG,AVZET,STDZET)
         CALL STATIS(CHI,NSEG,AVCHI,STDCHI)
         CALL STATIS(pseta,NSEG,AVpseta,STDpseta)
         CALL STATIS(pstheta,NSEG,AVpstheta,STDpstheta)
        IF(AVALF .GT. 180.0)AVALF = AVALF - 360.0
        IF(AVBET .GT. 180.0)AVBET = AVBET - 360.0
        IF(AVGAM .GT. 180.0)AVGAM = AVGAM - 360.0
        IF(AVDEL .GT. 180.0)AVDEL = AVDEL - 360.0
        IF(AVEPS .GT. 180.0)AVEPS = AVEPS - 360.0
        IF(AVZET .GT. 180.0)AVZET = AVZET - 360.0
        IF(AVCHI .GT. 180.0)AVCHI = AVCHI - 360.0
        IF(AVpseta .GT. 180.0)AVpseta = AVpseta - 360.0
        IF(AVpstheta .GT. 180.0)AVpstheta = AVpstheta - 360.0
       WRITE(11,120)AVALF,AVBET,AVGAM,AVDEL,AVEPS,AVZET,AVCHI,
     1AVpseta,AVpstheta
 120   FORMAT(/,7X,'Avg.',F7.1,8F9.1)
 121   FORMAT(7X,'S.D.',F7.1,8F9.1)
       WRITE(11,121)STDALF,STDBET,STDGAM,STDDEL,STDEPS,STDZET,STDCHI,
     1STDpseta,STDpstheta
C
        WRITE(11,114)
 114   FORMAT(//11X,' O4''-C1''  C1''-C2''  C2''-C3''  C3''-C4''  C4''-O
     14''    AMPL    PHASE   CLASS'/)
C
        SS3672 = SIN(36.0 * CONV) + SIN(72.0 * CONV)
C
         DO 2010 KS=1,NSEG
             TANPN = (ANEU4(KS) + ANEU1(KS)) - (ANEU3(KS) + ANEU0(KS))
           IF(ANEU2(KS).NE.0.000E0) THEN
             TANPD = 2.0 * ANEU2(KS) * SS3672
             PHASE(KS) = ATAN2(TANPN,TANPD)
             ANEUMAX(KS) = ANEU2(KS)/COS(PHASE(KS))
             PHASE(KS) = PHASE(KS) / CONV
           ELSE
             PHASE(KS) = 90.0
           END IF
C  --  PHASE ANGLE IS IN DEGREES NOW
         PH = PHASE(KS)
         IF(PH .LT. 0.0)PH = 360.0 + PHASE(KS)
         CLASS = EE(INT((PH + 9)/18) + 1)
         WRITE(11,113) ks,secnm(ks)(1:1),ANEU0(KS),ANEU1(KS),
     1     ANEU2(KS),ANEU3(KS),ANEU4(KS),ANEUMAX(KS),PHASE(KS),CLASS
 113    format('SG',i5,1X,a1,7F9.1,2X,A12)
 2010	CONTINUE
C
        CALL STATIS(ANEU0,NSEG,AVN0,STDN0)
        CALL STATIS(ANEU1,NSEG,AVN1,STDN1)
        CALL STATIS(ANEU2,NSEG,AVN2,STDN2)
        CALL STATIS(ANEU3,NSEG,AVN3,STDN3)
        CALL STATIS(ANEU4,NSEG,AVN4,STDN4)
        CALL STATIS(ANEUMAX,NSEG,AVNMAX,STDNMAX)
        CALL STATIS(PHASE,NSEG,AVPHS,STDPHS)
         WRITE(11,120)AVN0,AVN1,AVN2,AVN3,AVN4,AVNMAX,AVPHS
         WRITE(11,121)STDN0,STDN1,STDN2,STDN3,STDN4,STDNMAX,STDPHS
C
        RETURN
        END
      SUBROUTINE TORAN(A,B,C,D,CHI)
      DIMENSION A(3),B(3),C(3),D(3)
      DIMENSION WX(3),WY(3),WZ(3),EL(3),EM(3),EN(3)
C
      DO 100 J = 1,3
      WX(J)=A(J)-B(J)
      WY(J)=C(J)-B(J)
      WZ(J)=D(J)-B(J)
100   CONTINUE
      CALL UNITV(WX(1),WX(2),WX(3),WY(1),WY(2),WY(3),EL(1),EM(1),EN(1))
      CALL UNITV(WZ(1),WZ(2),WZ(3),WY(1),WY(2),WY(3),EL(2),EM(2),EN(2))
      CHI=ACOSF(EL(1)*EL(2)+EM(1)*EM(2)+EN(1)*EN(2))
      CALL UNITV(EL(2),EM(2),EN(2),EL(1),EM(1),EN(1),EL(3),EM(3),EN(3))
      CHECK=EL(3)*WY(1)+EM(3)*WY(2)+EN(3)*WY(3)
      IF(CHECK.LT.0.0) RETURN
      CHI=0.0-CHI
      IF(CHI.GT.180.0) THEN
          CHI = CHI - 360.0
      END IF
      RETURN
      END
C 
C***********************************************************************
C 
      SUBROUTINE UNITV(A1,A2,A3,B1,B2,B3,EL,EM,EN)
      C1 = A2*B3-A3*B2
      C2 = A3*B1-A1*B3
      C3 = A1*B2-A2*B1
      FAC = 0.0
      FAC1 =SQRTF(C1*C1+C2*C2+C3*C3)
      IF(FAC1.LE.0.00001) GO TO 100
      FAC = 1.0/FAC1
100   CONTINUE
      EL =C1*FAC
      EM =C2*FAC
      EN =C3*FAC
      RETURN
      END
C
      FUNCTION ACOSF(ANG)
1     FORMAT('  Cos (Theta) is greater than 1 !',F12.8)
      IF(ABS(ANG).LE.1.0) GO TO 102
C     WRITE(1,1) ANG
      IF (ABS(ANG).GT.1.000001) GO TO 102
      ACOSF=0.0
      IF (ANG.LT.0.0) ACOSF = 180.0
      RETURN
102   CONTINUE
      ACOSF=ACOS(ANG)*180.0/3.14159265
      RETURN
      END
C *********************************************************************
C
      FUNCTION SQRTF(X)
1     FORMAT(2X,'Square root of a negative quantity!',E16.8)
      IF (X.LT.0.0) GO TO 101
      SQRTF = SQRT(X)
      RETURN
101   CONTINUE
C     WRITE(1,1) X
      SQRTF = 0.0
      RETURN
      END
 
C***********************************************************************
C
      SUBROUTINE LINEFIT(DCS,NDATA,CENTRS)
	parameter ( nrs = 9999 )
	include 'parameters.h'
      DOUBLE PRECISION DCS,DCT
      DIMENSION DX(nrs),DY(nrs),DZ(nrs),AL(3),DIS(3),DCS(3),
     1          EL(3),EM(3),EN(3),CENTRS(nrs,3),DCT(3)
      CHARACTER*1 ANSO
C
C  Algorithm of this subroutine was adopted from K.-C. Chou, G. Nemethy
C  and H.A. Scheraga (1984) J. Amer. Chem. Soc. vol. 105, 3161-3170
C
      XSTAR = 0.0
      YSTAR = 0.0
      ZSTAR = 0.0
      IF(NDATA.GT.2) THEN
        DO 2000 I=1,NDATA
          XSTAR = XSTAR + CENTRS(I,1)
          YSTAR = YSTAR + CENTRS(I,2)
          ZSTAR = ZSTAR + CENTRS(I,3)
 2000	 CONTINUE
        XSTAR = XSTAR/NDATA
        YSTAR = YSTAR/NDATA
        ZSTAR = ZSTAR/NDATA
      AM11 = 0.0
      AM12 = 0.0
      AM13 = 0.0
      AM22 = 0.0
      AM23 = 0.0
      AM33 = 0.0
      DO 2001 I=1,NDATA
           DX(I) = CENTRS(I,1) - XSTAR
           DY(I) = CENTRS(I,2) - YSTAR
           DZ(I) = CENTRS(I,3) - ZSTAR
         AM11 = AM11 + DX(I) * DX(I)
         AM12 = AM12 + DX(I) * DY(I)
         AM13 = AM13 + DX(I) * DZ(I)
         AM22 = AM22 + DY(I) * DY(I)
         AM23 = AM23 + DY(I) * DZ(I)
         AM33 = AM33 + DZ(I) * DZ(I)
 2001	CONTINUE
      A = AM11
      B = AM12
      C = AM13
      D = AM22
      E = AM23
      F = AM33
      C1 = A*D*F + 2.0*B*C*E - A*E*E - B*B*F - C*C*D
      C2 = A*D + A*F + D*F - B*B - C*C - E*E
      C3 = A + D + F
      ALAMDA = 0.0
      ALH = 5.0
      DO 100 I=1,1000
      IF (ABS(ALH).GT.0.00001) THEN
      AFPR = (C2 + 2.0*ALAMDA*C3 - 3.0 * ALAMDA * ALAMDA)
      AF = C1 - ALAMDA*C2 + ALAMDA*ALAMDA*C3 - ALAMDA*ALAMDA*ALAMDA
      ALH = -AF/AFPR
      ALAMDA = ALAMDA + ALH
      ELSE
         GO TO 200
      END IF
 100  CONTINUE
 200  CONTINUE
C     WRITE(1,*) 'FIRST EIGEN-VALUE IS',ALAMDA
      A1 = 1.0
      A2 = ALAMDA - C3
      A3 = C1/ALAMDA
      SSQR = A2 * A2 - 4.0 * A3
      IF(SSQR.GT.0.0) THEN
         SSQR = SQRT(SSQR)
         ALAMD2 = (-A2 + SSQR) * 0.5
         ALAMD3 = (-A2 - SSQR) * 0.5
C        WRITE(1,*) 'SECOND EIGEN-VALUE IS',ALAMD2
C        WRITE(1,*) 'THIRD EIGEN-VALUE IS',ALAMD3
      ELSE
         SSQR = SQRT(-SSQR)
         A2 = -A2
C        WRITE(1,*) 'SECOND EIGEN-VALUE IS',A2,'i',SSQR
         SSQR = -SSQR
C        WRITE(1,*) 'THIRD EIGEN-VALUE IS',A2,'i',SSQR
      END IF
      AL(1) = ALAMDA
      AL(2) = ALAMD2
      AL(3) = ALAMD3
C     WRITE(1,*) 'EIGEN-VALUES',(AL(KKK),KKK=1,3)
      DO 2002 KT=1,3
         ADIS = 0.0
         A1 = AM11 - AL(KT)
         B1 = AM12
         C1 = AM13
         A2 = B1
         B2 = AM22 - AL(KT)
         C2 = AM23
C
         ELT = (B1 * C2 - C1 * B2)
         EMT = (C1 * A2 - A1 * C2)
         ENT = (A1 * B2 - B1 * A2)
         CONS = SQRT(ELT * ELT + EMT * EMT + ENT * ENT)
         ELT = ELT/CONS
         EMT = EMT/CONS
         ENT = ENT/CONS
C     WRITE(1,*) 'ELT,EMT,ENT',ELT,EMT,ENT,CONS
         DO 2003 K=1,NDATA
           ADIS = ADIS + DX(K)*DX(K) + DY(K)* DY(K) + DZ(K)*DZ(K)
           ADIS = ADIS - (ELT*DX(K)+EMT*DY(K)+ENT*DZ(K))**2
 2003	  CONTINUE
C     WRITE(1,*) 'MEAN DIST.',ADIS
         ADIS = SQRT(ADIS)
         DIS(KT) = ADIS
         EL(KT) = ELT
         EM(KT) = EMT
         EN(KT) = ENT
         IF(KT.GT.1) THEN
            IF(DIS(KT).LT.DIS(KT-1)) NN=KT
         END IF
 2002	CONTINUE
      DCS(1) = EL(NN) * 1.0
      DCS(2) = EM(NN) * 1.0
      DCS(3) = EN(NN) * 1.0
      ELSE
          DO 2004 KKN=1,3
            DCT(KKN) = CENTRS(1,KKN) - CENTRS(2,KKN)
 2004	   CONTINUE
          CALL NORMAL(DCT(1),DCT(2),DCT(3),DCS(1),DCS(2),DCS(3))
      END IF
      IF(CENTRS(1,3).GT.CENTRS(NDATA,3)) THEN
	   DCS(1) = -DCS(1)
	   DCS(2) = -DCS(2)
	   DCS(3) = -DCS(3)
      END IF
      RETURN
      END
C***********************************************************************
      FUNCTION ACOSD(X)
      IF(ABS(X).LT.1.0) THEN
         CONV=180.0/3.141592654
         ACOSD = ACOS(X) * CONV
      ELSE
C        WRITE(1,*) 'VALUE OF COS(THETA) > 1.0'
         ACOSD = 0.0
      END IF
      RETURN
      END
C
C-----------------------------------------------------------------------
C Subroutine to calc r, phi, z of P (and C1') atoms, dists etc.

        SUBROUTINE RTHEPH(ANSDB,ATTYPE,nb,ibapr)
	parameter ( natm = 999999 )
	parameter ( nrs = 9999 )
	include 'parameters.h'
	include 'coordinates.h'
      DOUBLE PRECISION XAXIS,YAXIS,ZAXIS,XM,YM,ZM1,ZM2,ZM3,Y1,Y2,Y3,
     1  Y12,Y22,Y32,X11,X12,X13,BMN,TL2,RL2,WTILT,WROLL,Y11,Y13,AMFO
        COMMON /BASINF/BASE(nrs),pair(3,nrs),prtype(nrs)
	DOUBLE PRECISION BPNORML,XH,YH,ZH,BPNORMLX,BPNORMLY,
     & MBPNORML,p1x,p1y,p1z,p2x,p2y,p2z,ax,ay,az,aux,auy,auz,
     & xc1n,yc1n,zc1n,xc1n1,yc1n1,zc1n1
	COMMON /BPZAXS/BPNORML(nrs,3),BPNORMLX(nrs,3),BPNORMLY(nrs,3),
     1     MBPNORML(nrs,3),AMFO(nrs,3)
	COMMON /GLOREO/DCRS1,DCRS2,DCRS3,ANGLE,DX1,DY1,INFBP(2,nrs)
      DIMENSION RAD(nrs),PHI(nrs),ZH(nrs),XH(nrs),YH(nrs),PPD1(nrs),
     1PPD2(50,nrs),NRES(nrs),DELPH(nrs),DELZ(nrs),ZP(nrs),Pc(nrs,3),
     2Pnu(3),YP(nrs),XP(nrs),ZPp(nrs),ibapr(nrs),xn(nrs),yn(nrs),
     3zn(nrs),xc1n(nrs),yc1n(nrs),zc1n(nrs),angc1n(nrs),
     4xc1n1(nrs),yc1n1(nrs),zc1n1(nrs)
      CHARACTER*4 attype,atnmlc
	character*1 ANSDB,ANS,ansnrm,pair,prtype,answw
      CHARACTER*3 BASE
      CONV=180.0/3.141592654
c
c  Copy relevent part of Coordinates from Common Array to local Phosphate array depending on
c  IBAPR information and leaving out the first Phosphate of each strand.
c
	IF(ATTYPE.EQ.'P   ') WRITE(11,1)
	IF(ATTYPE.EQ.'C1'' ') WRITE(11,2)
        L = 0
        L1=0
  	do ii=1,nb,2
	  iresimp=ibapr(ii)
	  nfoundtp=0
	  do j=ists(iresimp),iens(iresimp)
	    if(.not.(ii.eq.1.and.atnm(j).eq.'P   ')) then
	      if(atnm(j).eq.attype) then
	        L=L+1
                RAD(L) = SQRT(cr(j,1)*cr(j,1)+cr(j,2)*cr(j,2))
                PHI(L) = ATAN2(cr(j,2),cr(j,1)) * CONV
                ZH(L) = cr(j,3)
                XH(L) = cr(j,1)
                YH(L) = cr(j,2)
                NRES(L) = IRES(j)
	        nfoundtp=1
	      endif
c sukanya May 22, 2012
              if(attype.eq.'C1'' ')then
        if(((secnm(iresimp).eq.'ADE').or.(secnm(iresimp).eq.'GUA')).and.
     1  (atnm(j).eq.'N9 '))then
                   L1=L1+1
                   xn(L1)=cr(j,1)
                   yn(L1)=cr(j,2)
                   zn(L1)=cr(j,3)
        elseif(((secnm(iresimp).eq.'URA').or.(secnm(iresimp).eq.'CYT')
     1  .or.(secnm(iresimp).eq.'THY')).and.(atnm(j).eq.'N1 '))then
                   L1=L1+1
                   xn(L1)=cr(j,1)
                   yn(L1)=cr(j,2)
                   zn(L1)=cr(j,3)
        endif
              endif
	    endif
	  enddo
c	  if((ii.ne.1).and.(attype.eq.'C1'' '))then
	  if((attype.eq.'C1'' '))then
c            write(*,*)nres(l),xn(l),yn(l),zn(l),xh(l),yh(l),zh(l)
            xc1n1(L)=xn(l)-xh(l)
            yc1n1(L)=yn(l)-yh(l)
            zc1n1(L)=zn(l)-zh(l)
        call normal(xc1n1(l),yc1n1(l),zc1n1(l),xc1n(l),yc1n(l),zc1n(l))
          endif
	  if(nfoundtp.eq.0.and.ii.ne.1.and.iresimp.ne.0) then
	     write(6,*) 'An atom ',attype,' missing in Residue',iresimp
	     write(11,*) 'An atom ',attype,' missing in Residue',iresimp
	     l=l+1
	  endif
	enddo
	do i=nb,2,-2
	  iresimp=ibapr(i)
	  nfoundtp=0
	  do j=ists(iresimp),iens(iresimp)
	    if(.not.(i.eq.nb.and.atnm(j).eq.'P   ')) then
              if(atnm(j).eq.attype) then
                L=L+1
                RAD(L) = SQRT(cr(j,1)*cr(j,1)+cr(j,2)*cr(j,2))
                PHI(L) = ATAN2(cr(j,2),cr(j,1)) * CONV
                ZH(L) = cr(j,3)
                XH(L) = cr(j,1)
                YH(L) = cr(j,2)
                NRES(L) = IRES(j)
                nfoundtp=1
              endif
              if(attype.eq.'C1'' ')then
        if(((secnm(iresimp).eq.'ADE').or.(secnm(iresimp).eq.'GUA')).and.
     1  (atnm(j).eq.'N9 '))then
                   L1=L1+1
                   xn(L1)=cr(j,1)
                   yn(L1)=cr(j,2)
                   zn(L1)=cr(j,3)
        elseif(((secnm(iresimp).eq.'URA').or.(secnm(iresimp).eq.'CYT')
     1  .or.(secnm(iresimp).eq.'THY')).and.(atnm(j).eq.'N1 '))then
                   L1=L1+1
                   xn(L1)=cr(j,1)
                   yn(L1)=cr(j,2)
                   zn(L1)=cr(j,3)
        endif
              endif
	    endif
	  enddo
c	  if((i.ne.nb).and.(attype.eq.'C1'' ')) then
	  if((attype.eq.'C1'' ')) then
c            write(*,*)nres(l),xn(l),yn(l),zn(l),xh(l),yh(l),zh(l)
            xc1n1(L)=xn(l)-xh(l)
            yc1n1(L)=yn(l)-yh(l)
            zc1n1(L)=zn(l)-zh(l)
        call normal(xc1n1(l),yc1n1(l),zc1n1(l),xc1n(l),yc1n(l),zc1n(l))
          endif
	  if(nfoundtp.eq.0.and.i.ne.nb) then   !        .and.iresimp.ne.0) then
	     write(6,*) 'An atom ',attype,' missing in Residue',iresimp
	     write(11,*) 'An atom ',attype,' missing in Residue',iresimp
	     l=l+1
	  endif
	enddo
	Nphosp=L
c	write(*,*) natloc,' Atoms Now and',nb,' Important Residues in RTHEPH'
	write(9,*)'            B-P NORMAL(Z-Axis)','       P(x,y,z)I',
     &'              P(x,y,z)II','              PI - PII(unit vec)',
     &'         MEAN B-P ORIGIN','           Xp','    Yp','    Zp'

C
C    DELPHI, DELZ AND P--P DISTANCE CALC.
C 
         KK=Nphosp
         IF (ANSDB .EQ. 'Y') KK = Nphosp/2
         JJJ=1
         IF (ANSDB .EQ. 'Y') JJJ=2

	DO 2007 NP=1,Nphosp
	  IF(NP.LE.KK)THEN
C	vector(PI- PII) within a dinucleotide
	    ax = XH(NP) - XH(Nphosp - NP + 1)
	    ay = YH(NP) - YH(Nphosp - NP + 1)
	    az = ZH(NP) - ZH(Nphosp - NP + 1)
C	AMFO=avg. Middle Frame Base-pair origin
	    fp1x = XH(NP) - AMFO(NP,1)
	    fp1y = YH(NP) - AMFO(NP,2)
	    fp1z = ZH(NP) - AMFO(NP,3)
C	XP = [(PI - AMFO)*X] + [(PII - AMFO)*X] / 2.0
C	YP = [(PI - AMFO)*Y] + [(PII - AMFO)*-Y] / 2.0
C	ZP = [(PI - AMFO)*Z] + [(PII - AMFO)*-Z] / 2.0
	    fp2x = XH(Nphosp-NP+1) - AMFO(NP,1)
	    fp2y = YH(Nphosp-NP+1) - AMFO(NP,2)
	    fp2z = ZH(Nphosp-NP+1) - AMFO(NP,3)
C	
	 xp1=fp1x*BPNORMLX(NP,1)+fp1y*BPNORMLX(NP,2)+fp1z*BPNORMLX(NP,3)
	 xp2=fp2x*BPNORMLX(NP,1)+fp2y*BPNORMLX(NP,2)+fp2z*BPNORMLX(NP,3)
	    XP(NP) = 0.5*(xp1+xp2)

	    CALL NORMAL(ax,ay,az,aux,auy,auz)
            YP(NP)=0.5*(ax*BPNORMLY(NP,1)+ay*BPNORMLY(NP,2)
     &            +az*BPNORMLY(NP,3))
            ZP(NP)=0.5*(ax*BPNORML(NP,1)+ay*BPNORML(NP,2)
     &            +az*BPNORML(NP,3))
	
        write(9,3)NP,'.',BPNORML(NP,1),
     &  BPNORML(NP,2),BPNORML(NP,3),XH(NP),YH(NP),ZH(NP),XH(I-NP+1),
     &  YH(I-NP+1),ZH(I-NP+1),aux,auy,auz,AMFO(NP,1),AMFO(NP,2),
     &  AMFO(NP,3),XP(NP),YP(NP),ZP(NP)
C	write(*,*)NP,XP(NP),YP(NP),ZP(NP)
	  ENDIF
2007	CONTINUE

         I1 = 1
         I2 = KK - 1
C         I2 = KK
         IP = 0
C	write(*,*)KK, I2
       DO 2001 JJ = 1,JJJ
	WRITE(11,99) JJ
C          
	IF(ATTYPE.EQ.'P   ' .AND. JJ .EQ .2)
     1	WRITE(11,1001) NRES(IP1),RAD(IP1),PHI(IP1),ZH(IP1),
     2  XP(IP1),YP(IP1),ZP(IP1)

	IF(ATTYPE.EQ.'C1'' ' .AND. JJ .EQ .2)
     1	WRITE(11,981) NRES(IP1),RAD(IP1),PHI(IP1),ZH(IP1),
     2  XP(IP1),YP(IP1),ZP(IP1)
         DO 2002 NPP = I1,I2
          IP = IP + 1
	   NP = NPP
          NP1 = NP + 1
          IF(JJ .EQ. 2)NP1 = NP - 1
           DELP = PHI(NP1) - PHI(NP)
           IF(DELP .GT. 180.0)DELP = DELP - 360.0
           IF(DELP .LT. -180.0 .AND. DELP .GT. -360.0)DELP=360.0+DELP
           DELPH(NP) = DELP
           DELZ(NP) = ZH(NP1)  - ZH(NP)
           PPD1(NP)=SQRT((XH(NP)-XH(NP1))**2 + (YH(NP)-YH(NP1))**2 +
     1     (ZH(NP)-ZH(NP1))**2)
           ANGC1N(NP)=acos(xc1n(np)*xc1n(np1) + yc1n(np)*yc1n(np1) + 
     1     zc1n(np)*zc1n(np1))*conv
           

	IF(ATTYPE.EQ.'P   ') THEN
        WRITE(11,100)NRES(NP),RAD(NP),PHI(NP),ZH(NP),DELPH(NP),
     1  DELZ(NP),PPD1(NP),XP(NP),YP(NP),ZP(NP)

	ELSE

        WRITE(11,98) NRES(NP),RAD(NP),PHI(NP),ZH(NP),DELPH(NP),
     1  DELZ(NP),PPD1(NP),XP(NP),YP(NP),ZP(NP),angc1n(np)
	END IF
 2002	  CONTINUE
	IF(ATTYPE.EQ.'P   ' .AND. JJ .EQ. 1) 
     1  WRITE(11,1001) NRES(NP1),RAD(NP1),PHI(NP1),ZH(NP1),
     2  XP(NP1),YP(NP1),ZP(NP1)
C	IF(ATTYPE.EQ.'P   ' .AND. JJ .EQ. 1)write(*,*)NP1,ZP(NP1)
C	
	IF(ATTYPE.EQ.'C1'' '.AND. JJ .EQ. 1)
     1  WRITE(11,981) NRES(NP1),RAD(NP1),PHI(NP1),ZH(NP1),
     2  XP(NP1),YP(NP1),ZP(NP1)

C          I1 = KK + 1
C          I2 = Nphosp - 1
          I1 = KK + 2
	   I2 = Nphosp
          IP1= KK + 1
 2001	CONTINUE
C
         IK = IP
C	write(*,*)IK
        CALL STATIS(RAD,IK,AVRAD,STDRAD)
        CALL STATIS(DELPH,IK,AVPHI,STDPHI)
        CALL STATIS(DELZ,IK,AVZR,STDZR)
        CALL STATIS(PPD1,IK,AVPPD,STDPPD)
        CALL STATIS(XP,IK,XPAV,XPSTD)
        CALL STATIS(YP,IK,YPAV,YPSTD)
        CALL STATIS(ZP,IK,ZPAV,ZPSTD)
        call STATIS(ANGC1N,IK,AVC1N,STDC1N)
C	write(*,*)(ZP(NP),NP=1,LNDEX)
	 WRITE(11,111)AVRAD,AVPHI,AVZR,AVPPD,XPAV,
     2 YPAV,ZPAV,AVC1N
	 WRITE(11,112)STDRAD,STDPHI,STDZR,STDPPD,XPSTD,
     2 YPSTD,ZPSTD,STDC1N
C 
         IF (ANSDB .EQ. 'Y') THEN
C     CALCULATE INTERCHAIN P--P  and C1'--C1' DISTS.
C 
	IF(ATTYPE.EQ.'P   ') WRITE(11,101)
C---- C1'--C1' Interchain dist printout active.
	IF(ATTYPE.EQ.'C1'' ') WRITE(11,105)
C***********************first strnd #
        DO 2003 NP = 1,KK
         KI = KK + 1
         II=0
          DO 2004 NP2 = KI,Nphosp
           II = II + 1
	IF(ATTYPE.EQ.'P   ')THEN
           PPD2(NP,II)=SQRT((XH(NP)-XH(NP2))**2 + (YH(NP)-YH(NP2))**2 +
     1     (ZH(NP)-ZH(NP2))**2)
           PPDT(NP,II)=SQRT((XH(NP)-XH(NP2))**2 + (YH(NP)-YH(NP2))**2 +
     1     (ZH(NP)-ZH(NP2))**2)
		ENDIF
	IF(ATTYPE.EQ.'C1'' ')THEN
           PPD2(NP,II)=SQRT((XH(NP)-XH(NP2))**2 + (YH(NP)-YH(NP2))**2 +
     1     (ZH(NP)-ZH(NP2))**2)
           CCDT(NP,II)=SQRT((XH(NP)-XH(NP2))**2 + (YH(NP)-YH(NP2))**2 +
     1     (ZH(NP)-ZH(NP2))**2)
		ENDIF
C           PPD2(II)=SQRT((XH(NP)-XH(NP2))**2 + (YH(NP)-YH(NP2))**2 +
C     1     (ZH(NP)-ZH(NP2))**2)

 2004	  CONTINUE


C	K=0
C	J=0
	IF(ATTYPE.EQ.'P   ')THEN
	WRITE(11,102) NRES(NP),(PPDT(NP,IP),IP=1,II)
C	WRITE(11,102) NRES(NP),(PPD2(NP,IP),IP=1,II)
C	K=K+1
C	DO 1006 IP=1,II
C	J=J+1
C	PPDT(K,J)=PPD2(NP,IP)
C1006	CONTINUE
	ENDIF

	
C	M=0
C	L=0
	IF(ATTYPE.EQ.'C1'' ')THEN
	WRITE(11,104) NRES(NP),(CCDT(NP,IP),IP=1,II)
C	WRITE(*,104) NRES(NP),(PPD2(NP,IP),IP=1,II)
C	M=M+1
C        DO 1009 IP=1,II
C	L=L+1
C	CCDT(M,L)=PPD2(NP,IP)
CC	WRITE(*,'(11f8.2)')(CCDT(M,L)
C1009	CONTINUE
        ENDIF
C	DO 777 M=1,M
C	WRITE(*,'(11f8.2)')(CCDT(M,L),L=1,L)
C777	CONTINUE

 2003	CONTINUE
        ENDIF
         IF (ANSDB .EQ. 'Y'.AND.ATTYPE.EQ.'P   ') THEN
	DO 222 K=1,LNDEX
	 ZPp(K) = ZP(K)
222	CONTINUE
	ENDIF
c         WRITE(17,3312) FILENM,ZPAV
3312	FORMAT(1X,A20,f8.2)
C	IF(ANSDB .EQ. 'Y'.AND.ATTYPE.EQ.'P   ')WRITE(16,*)LNDEX
C         IF (ANSDB .EQ. 'Y') THEN
          IF (ANSDB .EQ. 'Y'.AND.ATTYPE.EQ.'C1'' ') THEN
C         IF (ANSDB .EQ. 'Y'.AND.ATTYPE.EQ.'P   ') THEN

	DO 2016 K=1,LNDEX

c         WRITE(12,212)BASE(K)(1:1),BASE(K+1)(1:1),INFBP(1,K),BASE(K),
c     1    INFBP(2,K),TILTL(K),ROLLL(K),TWISTL(K),SLX(K),SLY(K),SLZ(K),
c     2    APROP(K),DPROP(K),DZBUCK(K),DOPEN(K),ZPp(K)

c      if(K.eq.1)WRITE(16,126)BASE(K)(1:1),BASE(K+1)(1:1),INFBP(1,K),
c     1BASE(K),INFBP(2,K),TILTL(K),ROLLL(K),TWISTL(K),SLX(K),SLY(K),
c     2SLZ(K),PROTW(K),BUCAN(K),OPENAN(K),DZBUCK(K),ZPp(K),CCDT(K+2,II)
	
c      if(K.gt.1.and.K.lt.LNDEX-2)WRITE(16,126)
c     1BASE(K)(1:1),BASE(K+1)(1:1),INFBP(1,K),BASE(K),INFBP(2,K),
c     2TILTL(K),ROLLL(K),TWISTL(K),SLX(K),SLY(K),SLZ(K),PROTW(K),
c     3BUCAN(K),OPENAN(K),DZBUCK(K),ZPp(K),CCDT(K+2,II-K+1),
c     4PPDT(K+2,II-K+1),CCDT(K,II-K-2),PPDT(K-1,II-K-3)

c      if(K.eq.LNDEX-2)WRITE(16,126)
c     1BASE(K)(1:1),BASE(K+1)(1:1),INFBP(1,K),BASE(K),INFBP(2,K),
c     2TILTL(K),ROLLL(K),TWISTL(K),SLX(K),SLY(K),SLZ(K),PROTW(K),
c     3BUCAN(K),OPENAN(K),DZBUCK(K),ZPp(K),CCDT(K+2,II-K+1),
c     4PPDT(K+2,II-K+1),CCDT(K,II-K-2)

c      if(K.eq.LNDEX-1)WRITE(16,126)BASE(K)(1:1),BASE(K+1)(1:1),
c     1INFBP(1,K),BASE(K),INFBP(2,K),TILTL(K),ROLLL(K),TWISTL(K),
c     2SLX(K),SLY(K),SLZ(K),PROTW(K),BUCAN(K),OPENAN(K),DZBUCK(K),
c     3ZPp(K),CCDT(K,5)

c      if(K.eq.LNDEX)WRITE(16,126)BASE(K)(1:1),BASE(K+1)(1:1),
c     1INFBP(1,K),BASE(K),INFBP(2,K),TILTL(K),ROLLL(K),TWISTL(K),
c     2SLX(K),SLY(K),SLZ(K),PROTW(K),BUCAN(K),OPENAN(K),DZBUCK(K),
c     3ZPp(K)

 2016	CONTINUE
c      WRITE(12,200) INFBP(1,(LNDEX+1)),BASE(LNDEX+1),INFBP(2,(LNDEX+1))
c      WRITE(16,117) INFBP(1,(LNDEX+1)),BASE(LNDEX+1),INFBP(2,(LNDEX+1)),
c     2 PROTW(LNDEX+1),BUCAN(LNDEX+1),OPENAN(LNDEX+1)
	ENDIF
  
  1     FORMAT(//'  Phosphate polar coordinates & P--P distances'//
     116X,'R',6X,'Phi',8X,'Z',8X,'DelP',6X,'DelZ',8X,'P--P',6X,
     5 '   Xp','      Yp','       Zp',/)
  2   FORMAT(//'  C1'' cylindrical polar coordinates'//
     116X,'R',6X,'Phi',8X,'Z',8X,'DelP',6X,'DelZ',7X,'C1''--C1'' ',5X,
     1 'Xc1','     Yc1','     Zc1','  C1''-N9/N1'/)
3	format(i5,a1,3f9.4,6f8.3,6f9.4,3f6.2)
  98	FORMAT(1X,'C1''',I3,'.',3F9.2,2X,2F9.2,3X,F9.2,3X,4F9.2) 
 981	FORMAT(1X,'C1''',I3,'.',3F9.2,35X,3F9.2,11x) 
  99   FORMAT(3X,'Strand:',I2, ' 5''-->3'' ')
 100   FORMAT('ZP',1X,'  P',I3,'.',3F9.2,2X,2F9.2,3X,F9.2,3X,3F9.2)
1001   FORMAT('ZP',1X,'  P',I3,'.',3F9.2,35X,3F9.2)
 111   FORMAT(/4X,'AvCPC.',F7.2,20X,2F9.2,3X,F9.2,3X,4F9.2)
 112   FORMAT(4X,'S.D.',F9.2,20X,2F9.2,3X,F9.2,3X,4F9.2,/)
 126  FORMAT(A1,A1,1X,I2,2X,A3,2X,I2,15(1X,F7.2))
 117  FORMAT(3X,I2,2X,A3,2X,I2,49X,F7.2,1x,f7.2,1x,f7.2)
 101   FORMAT(/1X,' Inter-chain P--P distances'/'  Strand2: 5''-->3'' '
     1/' Strand1:')
 105   FORMAT(/1X,' Inter-chain C1''--C1'' distances'/'  Strand2: 5''-->
     13'' ' /' Strand1:')
 102   FORMAT('PP',2X,'P',I3,'.',20F7.2,2X,11(/,6X,20F7.2))
 104	FORMAT(2X,'C1''',I3,'.',20F7.2,11(/,9X,20F7.2))
 200  FORMAT(1X,I2,2X,A3,2X,I2,10(1X,F7.2))
 212  FORMAT(A1,A1,1X,I2,2X,A3,2X,I2,11(1X,F7.2))
        RETURN
        END
c
        SUBROUTINE FRACTION(A,B,C,ALPHA,BETA,GAMA)
	parameter ( natm = 999999 )
	parameter ( nrs = 9999 )
	include 'coordinates.h'
        DIMENSION AMAT(3,3),TEM(3)
	external matmul
C
        CONV=180.0/3.141592654
        ALPHA = ALPHA/CONV
        BETA = BETA/CONV        
        GAMA = GAMA/CONV
        AMAT(1,1)=1.0
        AMAT(1,2) = COS(GAMA)
        AMAT(1,3) = COS(BETA)
        AMAT(2,1) = 0.0
        AMAT(2,2) = SIN(GAMA)
        P1 = (COS(ALPHA) - COS(BETA) * COS(GAMA))/SIN(GAMA)     
        P2 = SQRT(SIN(BETA) * SIN(BETA) - P1 * P1)
        AMAT(2,3) = P1
        AMAT(3,1) = 0.0
        AMAT(3,2) = 0.0 
        AMAT(3,3) = P2
        DO 2000 KK=1,NAT
           CR(KK,1) = CR(KK,1) * A
           CR(KK,2) = CR(KK,2) * B
           CR(KK,3) = CR(KK,3) * C
           CALL MATMUL(AMAT,CR(KK,1),CR(KK,2),CR(KK,3),TEM)
           DO 2001 KKK=1,3
              CR(KK,KKK) = TEM(KKK)
 2001	    CONTINUE
 2000	CONTINUE
        RETURN
        END
C
      SUBROUTINE REORNT(IBAPR,CENTRS,NB,ANSORN,NPTS,ANSDB)
	parameter ( natm = 999999 )
	parameter ( nrs = 9999 )
	include 'coordinates.h'
      DIMENSION IBAPR(nrs),CENTRS(nrs,3),COOR1(3),COOR2(3)
      COMMON /DISPL/DXLOC(nrs),DYLOC(nrs),DZLOC(nrs),KKTI,
     1      BSCNTR(nrs,3),HORIG(nrs,3)
        CHARACTER*4 ATOM,ANSORN
        CHARACTER*1 ANSBC,ANSBB,ANSDB
C
 1      CONTINUE
        IF(ANSORN.EQ.'    ') THEN
         WRITE(8,91)
c 2016    WRITE(6,91)
  91    FORMAT(/' Should the GLOBAL AXIS be fitted to:'/' 1. Local HELIX
     1 ORIGINS (Type O)'/' 2. Base-Pair CENTERS (Type C)'/' 3. Any Backb
     2one atom (Type the atom name)')
c           READ(5,'(A4)') ATOM
           WRITE(8,'(5X,A4)') ansorn
c           ANSORN = ATOM
        END IF
C
        IF(ANSORN.EQ.'C   ') THEN
           DO 2000 KK=1,NB + 1
               DO 2001 KOO=1,3
                  CENTRS(KK,KOO) = BSCNTR(KK,KOO)
 2001	        CONTINUE
 2000	    CONTINUE
           NPTS = NB + 1
        ELSE IF(ANSORN.EQ.'O   ') THEN
           DO 2002 KS =1,NB
              DO 2003 KOO=1,3
                 CENTRS(KS,KOO)=HORIG(KS,KOO)
 2003	       CONTINUE
 2002	    CONTINUE
           NPTS = NB
        ELSE
          ATOM = ANSORN
	   ICNT = 0
          DO 2004 IS=1,NB
             IFIRST = IBAPR(IS*2-1)
             ISECND = IBAPR(IS*2)
             DO 2005 IATM=ISTS(IFIRST),IENS(IFIRST)
                IF(ATOM.EQ.ATNM(IATM)) THEN
	            ICNT = ICNT + 1
                   DO 2006 KOO=1,3
                      CENTRS(ICNT,KOO) = CR(IATM,KOO)
 2006	            CONTINUE
                END IF
 2005	      CONTINUE
	    IF(ISECND.NE.0) THEN
             DO 2007 IATM=ISTS(ISECND),IENS(ISECND)
                IF(ATOM.EQ.ATNM(IATM)) THEN
	            ICNT = ICNT + 1
                   DO 2008 KOO=1,3
                      CENTRS(ICNT,KOO) = CR(IATM,KOO)
 2008	            CONTINUE
                END IF
 2007	      CONTINUE
	    END IF
 2004	  CONTINUE
	  IF(ICNT.EQ.0) GO TO 2016
           ILAST = IBAPR((NB+1)*2-1)
           ILST = IBAPR((NB+1)*2)
           IF(IENS(ILST).LE.NAT.AND.ANSDB.EQ.'Y') THEN
              DO 2012 IATM=ISTS(ILAST),IENS(ILAST)
                IF(ATOM.EQ.ATNM(IATM)) THEN
	         ICNT = ICNT + 1
                   DO 2013 KOO=1,3
                      CENTRS(ICNT,KOO) = CR(IATM,KOO)
 2013	            CONTINUE
                END IF
 2012	      CONTINUE
	    IF(ILST.NE.0) THEN
             DO 2014 IATM=ISTS(ILST),IENS(ILST)
                IF(ATOM.EQ.ATNM(IATM)) THEN
	             ICNT = ICNT + 1
                   DO 2015 KOO=1,3
                      CENTRS(ICNT,KOO) = CR(IATM,KOO)
 2015	            CONTINUE
                END IF
 2014	      CONTINUE
	    END IF
          END IF
	   NPTS = ICNT
        END IF
        RETURN
2016	write(6,*) 'Invalid selection for molecular reorientation'
	stop
        END
C
	SUBROUTINE LINATM(DC,NDATA,CENTRS)
	parameter ( nrs = 9999 )
	DIMENSION CENTRS(nrs,3),VECTRS(3,nrs),DVA(nrs),NXA(nrs),DC(3)
	DOUBLE PRECISION DC,ELA,EMA,ENA
C
	DO 2000 KN=2,NDATA,2
	   NXA(KN) = KN
	   NXA(KN-1) = KN - 1
	   DO 2001 KOO=1,3
	      VECTRS(KOO,KN-1) = CENTRS(KN+1,KOO) - CENTRS(KN-1,KOO)
	      VECTRS(KOO,KN) = -CENTRS(KN,KOO) + CENTRS(KN+2,KOO)
 2001	   CONTINUE
 2000	CONTINUE
	if(ndata.ge.7) then
	  CALL PLANEo(VECTRS,(NDATA-2),NXA,0,DVA,ELA,EMA,ENA,DET,SDVY)
	else
	  write(6,801) (ndata-2)
	  stop
	endif
	DC(1) = ELA
	DC(2) = EMA
	DC(3) = ENA
	RETURN
 801	format('BAD Selection of BASE, only',I3,' atoms found')
	END
C
      SUBROUTINE PLANEo(X,NTOM,NX,NPLN,DV,EL,EM,EN,DET,SDV)
	parameter ( nrs = 9999 )
      DOUBLE PRECISION EL,EM,EN,Sl,Sm,Sn,y12(3),y13(3)
      DIMENSION X(3,nrs),NX(nrs),DV(nrs),Y(3,nrs),T(3,3),B(3),A(3)
      DIMENSION XY(3),XZ(3),VEC1(3),VEC2(3)
      external matmul
      INTEGER XY,XZ
C
        DMOVE = 0.0
	ntim = 0
C
 111   CONTINUE
        NCOUNT = 0
C
       DO 100 J1=1,NTOM
         NX1=NX(J1)
         IF(NX1.NE.0) THEN
            NCOUNT = NCOUNT + 1
            DO 2000 K=1,3
               Y(K,NCOUNT)=X(K,NX1) + DMOVE
 2000	     CONTINUE
         END IF
 100  CONTINUE
      if(ncount.ge.5) then
        NTOM1 = NCOUNT
      else
c	write(6,801) ncount
	stop 1
      endif
      DO 101 J1=1,3
      B(J1)=0.0
      DO 101 J2=1,3
      T(J1,J2)=0.0
 101  CONTINUE
      DO  102 J1=1,NTOM1
      T(1,1)=T(1,1)+Y(1,J1)*Y(1,J1)
      T(1,2)=T(1,2)+Y(1,J1)*Y(2,J1)
      T(1,3)=T(1,3)+Y(1,J1)*Y(3,J1)
      T(2,2)=T(2,2)+Y(2,J1)*Y(2,J1)
      T(2,3)=T(2,3)+Y(2,J1)*Y(3,J1)
      T(3,3)=T(3,3)+Y(3,J1)*Y(3,J1)
      B(1)=B(1)+Y(1,J1)
      B(2)=B(2)+Y(2,J1)
      B(3)=B(3)+Y(3,J1)
 102  CONTINUE
      DO 103 J1=1,3
      DO 103 J2=J1,3
      T(J2,J1)=T(J1,J2)
 103  CONTINUE
      CALL MINV(T,3,DX,XY,XZ)
      IF(DX.NE.0.0E0) THEN
      B1 = B(1)
      B2 = B(2)
      B3 = B(3)
      CALL MATMUL(T,B1,B2,B3,A)
      DET=1.0/(A(1)*A(1)+A(2)*A(2)+A(3)*A(3))
      DET=SQRT(DET)
      EL=A(1)*DET
      EM=A(2)*DET
      EN=A(3)*DET
      SDV = 0.0
      DO 104 J1=1,NTOM1
        DV(J1)=EL*Y(1,J1)+EM*Y(2,J1)+EN*Y(3,J1)-DET
   	 SDV = SDV + DV(J1) * DV(J1)    
 104  CONTINUE
      SDV = SQRT(SDV)
      IF(NPLN.NE.0.AND.SDV.GT.0.75) THEN
	   DMOVE = DMOVE + 25.0
	   ntim = ntim + 1
	   if(ntim.le.5) then
	     GO TO 111
	   else
C 
C LEAST SQUARE FIT did not work for the set of coordinates.  Going back to
C finding best plane through cross-product method and averaging.
C
	     sl=0.0D0
	     sm=0.0D0
	     sn=0.0D0
	     do 105 i=1,ntom1-2
	       do 106 k=1,3
	         y12(k) = y(k,i) - y(k,i+1)
		 y13(k) = y(k,i) - y(k,i+2)
 106	       continue
	       call cross(y12(1),y12(2),y12(3),y13(1),y13(2),y13(3),
     1  el,em,en)
     	       ang = (sl*el+sm*em+sn*en)/(sqrt(sl*sl+sm*sm+sn*sn)*
     1  sqrt(el*el+em*em+en*en))
     	       if(ang.lt.0.0) then
	         el = -el
		 em = -em
		 en = -en
		endif
		sl = sl+el
		sm = sm+em
		sn = sn+en
 105	      continue
	    endif
	 END IF
C
C  Commented out April 28, 2004 to run Unusual RNA structures also
C
c      ELSE
c      WRITE(6,*)' ERR ALL ATOMS ARE IN-PLANE, MOLECULE SHIFTED FOR CALC'
c         DMOVE = 25.0 + DMOVE
c         GO TO 111
      ENDIF
      RETURN
 801	format('BAD Selection of BASE, only',I3,' atoms found')
      END
C
C
	SUBROUTINE PONTSN(NDATA,CENTRS,DX,DY)
	parameter ( nrs = 9999 )
	DIMENSION CENTRS(nrs,3),POINT(2,nrs)
C
	DX = 0.0
	DY = 0.0
C
	DO 2000 I=1,(NDATA-2)
	  AM1=(CENTRS(I,2)-CENTRS(I+1,2))/(CENTRS(I,1)-CENTRS(I+1,1))
	  AM1 = -1.0/AM1
	  POI11 = (CENTRS(I,1) + CENTRS(I+1,1)) * 0.5
	  POI21 = (CENTRS(I,2) + CENTRS(I+1,2)) * 0.5
 	  AM2=(CENTRS(I+1,2)-CENTRS(I+2,2))/(CENTRS(I+1,1)-CENTRS(I+2,1))
	  AM2 = -1.0/AM2
	  POI12 = (CENTRS(I+1,1) + CENTRS(I+2,1)) * 0.5
	  POI22 = (CENTRS(I+1,2) + CENTRS(I+2,2)) * 0.5
	  C1 = POI21 - AM1 * POI11
	  C2 = POI22 - AM2 * POI12
	  POINT(1,I) = (C1 - C2)/(AM2 - AM1)
	  POINT(2,I) = (AM2 * C1 - AM1 * C2)/(AM2 - AM1)
	  DX = DX + POINT(1,I)
	  DY = DY + POINT(2,I)
 2000	CONTINUE
	DX = DX/(NDATA-2)
	DY = DY/(NDATA-2)
	RETURN
	END
C
	SUBROUTINE PONTDB(NDATA,CENTRS,DX,DY)
	parameter ( nrs = 9999 )
	DIMENSION CENTRS(nrs,3),POINT(2,nrs)
C
	DX = 0.0
	DY = 0.0
	NVEC = NDATA/2
C
	DO 2000 I=1,(NDATA-2),2
	  AM1=(CENTRS(I,2)-CENTRS(I+1,2))/(CENTRS(I,1)-CENTRS(I+1,1))
	  AM1 = -1.0/AM1
	  POI11 = (CENTRS(I,1) + CENTRS(I+1,1)) * 0.5
	  POI21 = (CENTRS(I,2) + CENTRS(I+1,2)) * 0.5
 	  AM2=(CENTRS(I+2,2)-CENTRS(I+3,2))/(CENTRS(I+2,1)-CENTRS(I+3,1))
	  AM2 = -1.0/AM2
	  POI12 = (CENTRS(I+2,1) + CENTRS(I+3,1)) * 0.5
	  POI22 = (CENTRS(I+2,2) + CENTRS(I+3,2)) * 0.5
	  C1 = POI21 - AM1 * POI11
	  C2 = POI22 - AM2 * POI12
	  POINT(1,I) = (C1 - C2)/(AM2 - AM1)
	  POINT(2,I) = (AM2 * C1 - AM1 * C2)/(AM2 - AM1)
	  DX = DX + POINT(1,I)
	  DY = DY + POINT(2,I)
 2000	CONTINUE
	DX = DX/(NVEC - 1)
	DY = DY/(NVEC - 1)
	RETURN
	END
C
	SUBROUTINE LINASN(DC,NDATA,CENTRS)
	parameter ( nrs = 9999 )
	DIMENSION CENTRS(nrs,3),VECTRS(3,nrs),DVA(nrs),NXA(nrs),DC(3)
	DOUBLE PRECISION DC,ELA,EMA,ENA
C
	DO 2000 KN=2,NDATA
	   NXA(KN-1) = KN - 1
	   DO 2001 KOO=1,3
	      VECTRS(KOO,KN-1) = CENTRS(KN,KOO) - CENTRS(KN-1,KOO)
 2001	   CONTINUE
 2000	CONTINUE
	if(ndata.ge.6) then
	  CALL PLANEo(VECTRS,(NDATA-1),NXA,0,DVA,ELA,EMA,ENA,DET,SDVZ)
	else
c	  write(6,801) (ndata-1)
	  stop
	endif
	DC(1) = ELA
	DC(2) = EMA
	DC(3) = ENA
	RETURN
 801	format('BAD Selection of BASE, only',I3,' atoms found')
	END

	subroutine findybs(ibase,xaxis,yaxis,zaxis,angle,yo,
     1                     type,cystrns)
	parameter ( natm = 999999 )
	parameter ( nrs = 9999 )
	include 'coordinates.h'
	dimension yaxis(3),yo(3),yp(3),puratm(9),pyratm(6),xyzbs(3,9),
     1   nbase(9),dva(9),nbs(9),xaxis(3),zaxis(3),xch(3)
	double precision el,em,en,xaxis,yaxis,zaxis
	character*4 puratm,pyratm
	character*1 type,cystrns
	logical purine
C
	data puratm/'N9  ','C8  ','N7  ','C6  ','C5  ','C4  ','N3  ',
     1  'C2  ','N1  '/
     	data pyratm/'N1  ','C6  ','C5  ','C4  ','N3  ','C2  '/
	do i=1,9
	  nbs(i)=0
	enddo
	purine=.false.
c	write(*,'(i5,2x,a3,2x,a1)') ibase,secnm(ibase),type
	do while(secnm(ibase)(1:1).eq.' ') 
	  secnm(ibase)(1:1)=secnm(ibase)(2:2)
	  secnm(ibase)(2:2)=secnm(ibase)(3:3)
	  secnm(ibase)(3:3)= ' '
	enddo

c	write(*,*) ibase,ists(ibase),iens(ibase),secnm(ibase)
c	write(9,1) secnm(ibase),ibase,type,cystrns
1	format('Res.',a3,' No.',i5,' pairing with ',a1,' face in ',a1,
     1 ' orientation')
	  if(secnm(ibase).eq.'THY'.or.secnm(ibase).eq.'CYT'.or.
     1  secnm(ibase).eq.'T  '.or.secnm(ibase).eq.'C  '.or.secnm(ibase)
     2  .eq.'U  '.or.secnm(ibase).eq.'URA') then
c     	write(9,*) 'Pyrimidine found',ibase,secnm(ibase)
	    nsel = 0
	    do iatm=ists(ibase),iens(ibase)
	     if(type.ne.'H'.and.type.ne.'h'.and.type.ne.'S'.and.
     1          type.ne.'s'.and.type.ne.'z') then
	      if(atnm(iatm).eq.'N3  ') then
	         yo(1)=cr(iatm,1)
		 yo(2)=cr(iatm,2)
		 yo(3)=cr(iatm,3)
	      elseif(atnm(iatm).eq.'O4  ') then
	         yp(1)=cr(iatm,1)
		 yp(2)=cr(iatm,2)
		 yp(3)=cr(iatm,3)
	      elseif(atnm(iatm).eq.'N4  ') then
	         yp(1)=cr(iatm,1)
		 yp(2)=cr(iatm,2)
		 yp(3)=cr(iatm,3)
	      elseif(atnm(iatm).eq.'C6  ') then
	         xch(1)=cr(iatm,1)
	         xch(2)=cr(iatm,2)
	         xch(3)=cr(iatm,3)
	      endif
	     elseif(type.eq.'H'.or.type.eq.'h') then
	      if(atnm(iatm).eq.'C5  ') then
	        yo(1)=cr(iatm,1)
	        yo(2)=cr(iatm,2)
	        yo(3)=cr(iatm,3)
	      elseif(atnm(iatm).eq.'O4  '.or.atnm(iatm).eq.'N4  ') then
	        yp(1)=cr(iatm,1)
	        yp(2)=cr(iatm,2)
	        yp(3)=cr(iatm,3)
	      elseif(atnm(iatm).eq.'C2  ') then
	        xch(1)=cr(iatm,1)
	        xch(2)=cr(iatm,2)
	        xch(3)=cr(iatm,3)
	      endif
	     elseif(type.eq.'S'.or.type.eq.'z'.or.type.eq.'s') then
	      if(atnm(iatm).eq.'C1* '.or.atnm(iatm).eq.'C1'' ') then
	        yo(1)=cr(iatm,1)
	        yo(2)=cr(iatm,2)
	        yo(3)=cr(iatm,3)
	      elseif(atnm(iatm).eq.'O2  ') then
	        yp(1)=cr(iatm,1)
	        yp(2)=cr(iatm,2)
	        yp(3)=cr(iatm,3)
	      elseif(atnm(iatm).eq.'C6  ') then
	        xch(1)=cr(iatm,1)
	        xch(2)=cr(iatm,2)
	        xch(3)=cr(iatm,3)
	      endif
	     endif

	      do ipy=1,6
	        if(atnm(iatm).eq.pyratm(ipy).and.nsel.lt.6) then
		  nsel=nsel+1
		  xyzbs(1,nsel)=cr(iatm,1)
		  xyzbs(2,nsel)=cr(iatm,2)
		  xyzbs(3,nsel)=cr(iatm,3)
		  nbs(nsel)=nsel
		endif
	      enddo
	    enddo
c	      write(*,*) 'Selected atoms for pyrimidine',nsel

c	    write(6,*) 'Pyr ',ibase,secnm(ibase),(yo(kx),kx=1,3),(yp(kx),
c     1  kx=1,3)
	    xch(1)=yo(1)-xch(1)
	    xch(2)=yo(2)-xch(2)
	    xch(3)=yo(3)-xch(3)
	  elseif(secnm(ibase).eq.'ADE'.or.secnm(ibase).eq.'GUA'
     1   .or.secnm(ibase).eq.'A  '.or.secnm(ibase).eq.'G  '.or.
     2   secnm(ibase).eq.'I  ') then
c     	write(9,*) 'Purine found',ibase,secnm(ibase)
	    purine=.true.
	    nsel = 0
	    do iatm=ists(ibase),iens(ibase)
	      if(type.ne.'H'.and.type.ne.'h'.and.type.ne.'S'.and.
     1    type.ne.'s'.and.type.ne.'z') then
	        if(atnm(iatm).eq.'N1  ') then
	          yo(1)=cr(iatm,1)
	          yo(2)=cr(iatm,2)
	          yo(3)=cr(iatm,3)
	        elseif(atnm(iatm).eq.'N6  ') then
	          yp(1)=cr(iatm,1)
	          yp(2)=cr(iatm,2)
	          yp(3)=cr(iatm,3)
	        elseif(atnm(iatm).eq.'O6  ') then
	          yp(1)=cr(iatm,1)
	          yp(2)=cr(iatm,2)
	          yp(3)=cr(iatm,3)
	        elseif(atnm(iatm).eq.'N9  ') then
	          xch(1)=cr(iatm,1)
	          xch(2)=cr(iatm,2)
	          xch(3)=cr(iatm,3)
	        endif
	      elseif(type.eq.'H'.or.type.eq.'h') then
	        if(atnm(iatm).eq.'N7  ') then
	          yo(1)=cr(iatm,1)
	          yo(2)=cr(iatm,2)
	          yo(3)=cr(iatm,3)
	        elseif(atnm(iatm).eq.'O6  '.or.atnm(iatm).eq.'N6  ') then
	          yp(1)=cr(iatm,1)
	          yp(2)=cr(iatm,2)
	          yp(3)=cr(iatm,3)
	        elseif(atnm(iatm).eq.'C4  ') then
	          xch(1)=cr(iatm,1)
	          xch(2)=cr(iatm,2)
	          xch(3)=cr(iatm,3)
	        endif
	      elseif(type.eq.'S'.or.type.eq.'s'.or.type.eq.'z') then
	        if(atnm(iatm).eq.'C1* '.or.atnm(iatm).eq.'C1'' ') then
	          yp(1)=cr(iatm,1)
	          yp(2)=cr(iatm,2)
	          yp(3)=cr(iatm,3)
	        elseif(atnm(iatm).eq.'N3  ') then
	          yo(1)=cr(iatm,1)
	          yo(2)=cr(iatm,2)
	          yo(3)=cr(iatm,3)
	        elseif(atnm(iatm).eq.'C5  ') then
	          xch(1)=cr(iatm,1)
	          xch(2)=cr(iatm,2)
	          xch(3)=cr(iatm,3)
	        endif
              endif 
	      do ipu=1,9
	        if(atnm(iatm).eq.puratm(ipu).and.nsel.lt.9) then
		  nsel=nsel+1
		  xyzbs(1,nsel)=cr(iatm,1)
		  xyzbs(2,nsel)=cr(iatm,2)
		  xyzbs(3,nsel)=cr(iatm,3)
		  nbs(nsel)=nsel
		endif
	      enddo
	    enddo
c              write(*,*) 'Selected atoms for purine',nsel
c	    write(6,*) 'Pur ',ibase,secnm(ibase),(yo(kx),kx=1,3),(yp(kx)
c     1  ,kx=1,3)
	    xch(1)=yo(1)-xch(1)
	    xch(2)=yo(2)-xch(2)
	    xch(3)=yo(3)-xch(3)
	  endif
	  alen=0.0
	  do kx=1,3
	    yaxis(kx)=yp(kx)-yo(kx)
	    alen=alen+yaxis(kx)*yaxis(kx)
	  enddo
	  if(cystrns.eq.'T')then
	    do kx=1,3
	       yaxis(kx)= -yaxis(kx)
	    enddo
	  endif
	  if(purine.and.(type.eq.'S'.or.type.eq.'z'.or.type.eq.'s')) then
	    do kx=1,3
	       yaxis(kx)= -yaxis(kx)
	    enddo
	  endif
	      
	  alen = sqrtf(alen)
	  do kx=1,3
	    yaxis(kx)=yaxis(kx)/alen
	  enddo
	  if(nsel.ge.6) then
c	    write(*,*) nbs
	    call planeb(xyzbs,nsel,nbs,1,dva,el,em,en,detr,sdv1)
c	    write(*,*) el,em,en
	    xaxis(1)=el
	    xaxis(2)=em
	    xaxis(3)=en
C 
C  Base X-axis is now the vector perpendicular to the plane of the BASE
C
	    call cross(el,em,en,yaxis(1),yaxis(2),yaxis(3),
     1 zaxis(1),zaxis(2),zaxis(3))
            angle = zaxis(1)*xch(1)+zaxis(2)*xch(2)+zaxis(3)*xch(3)
	    alen=sqrt(xch(1)*xch(1)+xch(2)*xch(2)+xch(3)*xch(3))
	    angle=angle/alen
	    if(angle.gt.0.0) then
	       zaxis(1)=-zaxis(1)
	       zaxis(2)=-zaxis(2)
	       zaxis(3)=-zaxis(3)
	       xaxis(1)=-xaxis(1)
	       xaxis(2)=-xaxis(2)
	       xaxis(3)=-xaxis(3)
	    endif
	  endif 

	write(9,*) 'Base No.',ibase,' Xaxis=',xaxis,' Yaxis=',yaxis,
     1             ' Zaxis=',zaxis
987     format('Base No.  ',i3,3(a7,3f10.6))
	write(21,987) ibase,' Xaxis=',xaxis,' Yaxis=',yaxis,
     1             'Zaxis=',zaxis
c	write(6,*) 'Base No.',ibase,' Xaxis=',xaxis,' Yaxis=',yaxis,
c     1             ' Zaxis=',zaxis
	  return
	  end

C --------------------------------------------------------------------
C    SUBROUTINE FOR TIP AND INCL. CALC. W.R.T. LOCAL HELIX AXIS &
C     WTILT AND WROLL CALC. W.R.T. THE MEAN DOUBLET AXES.
C --------------------------------------------------------------------
      SUBROUTINE WEGDMbs(XAXIS,YAXIS,BMN,param,STKINV)
	parameter ( nrs = 9999 )
	include 'parameters.h'
      DOUBLE PRECISION XAXIS,YAXIS,ZAXIS,XM,YM,ZM1,ZM2,ZM3,Y1,Y2,Y3,
     1   Y12,Y22,Y32,X11,X12,X13,TL2,RL2,WTILT,WROLL,Y11,Y13
      LOGICAL STKINV,RETCAL
      DIMENSION XAXIS(2,3),YAXIS(2,3),XM(3),YM(3),ZAXIS(3),BMN(3)
     1       ,param(6)
C
      CONV=180.0/3.141592
C
      XM(1) = XAXIS(1,1) + XAXIS(2,1)
      XM(2) = XAXIS(1,2) + XAXIS(2,2)
      XM(3) = XAXIS(1,3) + XAXIS(2,3)
      YM(1) = YAXIS(1,1) + YAXIS(2,1)
      YM(2) = YAXIS(1,2) + YAXIS(2,2)
      YM(3) = YAXIS(1,3) + YAXIS(2,3)
c      write(9,*) 'Xm=',xm
c      write(9,*) 'Ym=',ym
      CALL CROSS(XM(1),XM(2),XM(3),YM(1),YM(2),YM(3),ZM1,ZM2,ZM3)
      CALL NORMAL(ZM1,ZM2,ZM3,ZAXIS(1),ZAXIS(2),ZAXIS(3))
      CALL NORMAL(XM(1),XM(2),XM(3),X11,X12,X13)
      CALL NORMAL(YM(1),YM(2),YM(3),Y11,Y12,Y13)
C
      XM(1) = X11
      XM(2) = X12
      XM(3) = X13
      YM(1) = Y11
      YM(2) = Y12
      YM(3) = Y13
      write(9,1) xm
1	format(' MEAN X-AXIS for basepair',1X,F12.8,1x,F12.8,1X,F12.8)
      write(9,2) ym
2	format(' MEAN Y-AXIS for basepair',1X,F12.8,1x,F12.8,1X,F12.8)
      angxmym=acos(x11*y11+x12*y12+x13*y13)*conv
      write(9,*) 'Angle between Xm & Ym= for basepair',angxmym
      TL2 = 0.0
      RL2 = 0.0
      CNEW = 0.0
      TSLIDE = 0.0
      SLIDE = 0.0
      WRITE(9,3) (ZAXIS(I),I=1,3)
3	format(' MEAN Z-AXIS for basepair',1X,F12.8,1x,F12.8,1X,F12.8)
      DO 2002 I=1,3
         TL2 = TL2 - YAXIS(1,I) * ZAXIS(I)
         RL2 = RL2 + XAXIS(1,I) * ZAXIS(I)
         CNEW = CNEW + XM(I) * BMN(I)
         TSLIDE = TSLIDE + BMN(I) * ZAXIS(I)
         SLIDE = SLIDE + BMN(I) * YM(I)
 2002 CONTINUE
C
      WTILT =  2.0 * ASIN(TL2) * CONV
      WROLL =  2.0 * ASIN(RL2) * CONV
C
C Signs of BUCKLE (param(2)) and STAGGER (param(5) are inverted to make it 
C consistent with X3DNA definitions. 		bhatta, June 18, 2004
C
      param(1)=wtilt
      param(2)=-wroll
      param(4)=cnew
      param(5)=-slide
      param(6)=tslide
      CALL CROSS(XAXIS(1,1),XAXIS(1,2),XAXIS(1,3),ZAXIS(1),ZAXIS(2),ZAXI
     1S(3),Y1,Y2,Y3)
      CALL CROSS(XAXIS(2,1),XAXIS(2,2),XAXIS(2,3),ZAXIS(1),ZAXIS(2),ZAXI
     1S(3),Y12,Y22,Y32)
      ANG = Y1 * Y12 + Y2 * Y22 + Y3 * Y32
      ANGLE = ACOS(ANG) * CONV
C
      CALL CROSS(Y1,Y2,Y3,Y12,Y22,Y32,ZM1,ZM2,ZM3)
      ANG = ZM1 * ZAXIS(1) + ZM2 * ZAXIS(2) + ZM3 * ZAXIS(3)
      IF(ANG.LT.0.0) ANGLE = - ANGLE
C
c	TWISTL(LNDEX) = ANGLE
	param(3)=angle
c      CALL SLIDCAL
c      IF(STKINV) THEN
c  	  DO 3000 KK=1,3
c	    XAXIS(2,KK) = -XAXIS(2,KK)
c	    YAXIS(2,KK) = -YAXIS(2,KK)
c 3000	  CONTINUE
c      END IF
C
C
c	write(*,*) (param(k6),k6=1,6)
      RETURN
      END
C
      SUBROUTINE PLANEB(X,NTOM,NX,NPLN,DV,EL,EM,EN,DET,SDV)
	parameter ( nrs = 9999 )
      COMMON /REPLY/ANSOR,ANSC1,ANSBPN,ANSDB,ANSrnt,ANSHO,ansnrm,answw,
     1   anspp,anstor,ansovl,anshlx
      DOUBLE PRECISION EL,EM,EN,Sl,Sm,Sn,y12,y13
      DIMENSION X(3,nrs),NX(nrs),DV(nrs),Y(3,nrs),T(3,3),B(3),A(3)
      DIMENSION XY(3),XZ(3),VEC1(3),VEC2(3),y12(3),y13(3)
      external matmul
      CHARACTER*1 ANSOR,ANSC1,ANSrnt,ANSHO,ANSDB,ANSBPN,pair,ansnrm,
     1        answw,anspp,anstor,ansovl,anshlx
      INTEGER XY,XZ
C
      DMOVE = 0.0
      ntim = 0
      NCOUNT = 0
C      write(*,*)NTOM
      DO 100 J1=1,NTOM
        NX1=NX(J1)
        IF(NX1.NE.0) THEN
          NCOUNT = NCOUNT + 1
          DO 2000 K=1,3
C      write(*,*)J1,NX1,NCOUNT,'  Working-Sukanya'
            Y(K,NCOUNT)=X(K,NX1) 
2000      CONTINUE
        END IF
100   CONTINUE

      if(ncount.ge.4) then
        NTOM1 = NCOUNT
      else
	write(40,*) ncount,ntom,nx1,' nx=',nx
c	write(6,801) ncount
	stop 1
      endif

       if(ansnrm.ne.'Y') then
 	 do while(ntim.le.5)
           DO 101 J1=1,3
             B(J1)=0.0
             DO 101 J2=1,3
             T(J1,J2)=0.0
 101       CONTINUE
           do j1=1,ntom1
             do k=1,3
               y(k,j1)=y(k,j1)+0.003*(rand(0)-0.5)
             enddo
c             write(*,*) y(3,j1)
           enddo
           DO  102 J1=1,NTOM1
             T(1,1)=T(1,1)+Y(1,J1)*Y(1,J1)
             T(1,2)=T(1,2)+Y(1,J1)*Y(2,J1)
             T(1,3)=T(1,3)+Y(1,J1)*Y(3,J1)
             T(2,2)=T(2,2)+Y(2,J1)*Y(2,J1)
             T(2,3)=T(2,3)+Y(2,J1)*Y(3,J1)
             T(3,3)=T(3,3)+Y(3,J1)*Y(3,J1)
             B(1)=B(1)+Y(1,J1)
             B(2)=B(2)+Y(2,J1)
             B(3)=B(3)+Y(3,J1)
 102       CONTINUE
           DO 103 J1=1,3
             DO 103 J2=J1,3
               T(J2,J1)=T(J1,J2)
 103       CONTINUE
           CALL MINV(T,3,DX,XY,XZ)
	   write(*,*) dx
           IF(DX.NE.0.0E0) THEN
             B1 = B(1)
             B2 = B(2)
             B3 = B(3)
             CALL MATMUL(T,B1,B2,B3,A)
             DET=1.0/(A(1)*A(1)+A(2)*A(2)+A(3)*A(3))
             DET=SQRT(DET)
             EL=A(1)*DET
             EM=A(2)*DET
             EN=A(3)*DET
             SDV = 0.0
             DO 104 J1=1,NTOM1
               DV(J1)=EL*Y(1,J1)+EM*Y(2,J1)+EN*Y(3,J1)-DET
   	       SDV = SDV + DV(J1) * DV(J1)    
 104         CONTINUE
             SDV = SQRT(SDV)
             IF(NPLN.NE.0.AND.SDV.GT.0.75) THEN
c	       DMOVE = DMOVE + 25.0
	       ntim = ntim + 1
c	       write(*,*) 'Ntime=',ntim, npln, sdv
	       do j1=1,ntom1
	         do k=1,3
	           y(k,j1)=y(k,j1)+0.1*(rand(0)-0.5)
	         enddo
c		 write(*,*) y(3,j1)
	       enddo
	     else
	       return
	     endif
	   endif
	 enddo

	else
C 
C LEAST SQUARE FIT did not work for the set of coordinates.  Going back to
C finding best plane through cross-product method and averaging.
C
c           write(*,*) 'Cross-product method with',ntom1-2
	   do k=1,3
	      y12(k)=y(k,1) - y(k,2)
	      y13(k)=y(k,1) - y(k,3)
	   enddo

	   call cross(y12(1),y12(2),y12(3),y13(1),y13(2),y13(3),
     1 sl,sm,sn)
c       write(*,*) 'sl',sl,sm,sn
        do i=1,ntom1
          do j=i+1,ntom1
            do k=1,ntom1
              if(i.ne.j.and.i.ne.k.and.j.ne.k) then
                if(j.eq.2.and.i.eq.1.and.k.eq.3) then
                  xx=1
                elseif(j.eq.2.and.i.eq.3.and.k.eq.1) then
                  xx=1
		else
c		write(*,*) i,j,k
	        do 106 ll=1,3
	          y12(ll) = y(ll,i) - y(ll,j)
	          y13(ll) = y(ll,k) - y(ll,j)
 106	        continue
                call cross(y12(1),y12(2),y12(3),y13(1),y13(2),y13(3),
     1  el,em,en)
c                write(*,*) i,j,k,el,em,en
     	        ang = (sl*el+sm*em+sn*en)/(sqrt(sl*sl+sm*sm+sn*sn)*
     1  sqrt(el*el+em*em+en*en))
     	        IF(ang.lt.0.0) then
	           el = -el
	           em = -em
	           en = -en
	        endif
  	        sl = sl+el
	        sm = sm+em
	        sn = sn+en
c	       write(*,*) 'Added sl',sl,sm,sn
               endif
	      endif
	    enddo
	   enddo
	  enddo
          call normal(sl,sm,sn,el,em,en)
c	  write(*,*) el,em,en, 'at the end'

	  endif
      RETURN
 801	format('BAD Selection of BASE, only',I3,' atoms found in PLANEB')
      END
	subroutine getblprm(ibp1,ibp2,typeadd,bptadd,
     2 bptilt,bproll,bptwst,bpshft,bpslid,bprise,stk,ybs)
	parameter ( natm = 999999 )
	parameter ( nrs = 9999 )
	include 'coordinates.h'
	dimension xbase1(3),ybase1(3),zbase1(3),xbase2(3),ybase2(3),
     1            zbase2(3),orig1(3),orig2(3),bm(3),xbs(2,3),ybs(2,3),
     2            param(6)
	double precision xbase1,ybase1,zbase1,xbase2,ybase2,zbase2,xbs,
     1            ybs
	character*3 typeadd
	character*1 bptadd,type1,type2
	logical stk
	
c	write(*,*) 'In GETBLPRM',ibp1,ibp2
	type1 = typeadd(1:1)
	type2 = typeadd(3:3)
	
c	write(*,*) ibp1,ibp2,type1,type2,bptadd,secnm(ibp1),secnm(ibp2)
          call findybs(ibp1,xbase1,ybase1,zbase1,ang1,
     1                 orig1,type1,bptadd)
c	write(*,*) 'After first call to findybs',ibp1,ibp2,type1,type2
          call findybs(ibp2,xbase2,ybase2,zbase2,ang2,
     1                 orig2,type2,'C')
c	write(*,*) 'After 2nd call to findybs',ibp1,ibp2,type1,type2
	  ang=0.0
          do kx=1,3
            ang=ang+zbase1(kx)*zbase2(kx)
            bm(kx)=orig2(kx)-orig1(kx)
            xbs(1,kx)=-xbase1(kx)
            xbs(2,kx)=xbase2(kx)
            ybs(1,kx)=ybase1(kx)
            ybs(2,kx)=ybase2(kx)
          enddo
          ang=acosf(ang)
c         write(9,'(3f10.6)') ((xbs(kn,kx),kx=1,3),kn=1,2)
c         write(9,'(3f10.6)') ((ybs(kn,kx),kx=1,3),kn=1,2)
c	 write(9,'(3f10.6)') (-zbase1(kx),kx=1,3)
c	 write(9,'(3f10.6)') (zbase2(kx),kx=1,3)
c         write(*,*) 'bm=',bm
	  stk=.false.
	  write(9,1) orig1
 	  write(9,1) orig2
!1	format('ORIGIN',3f8.2)						!modified by DM
	  write(21,1) orig1
 	  write(21,1) orig2 	  
1	format('Base ORIGIN',3f8.3)
	write(9,*)'Xbs & Ybs1 b4 WEGDMbs',(xbs(1,Lc),ybs(1,Lc),Lc=1,3)
	write(9,*)'Xbs & Ybs2 b4 WEGDMbs',(xbs(2,Lc),ybs(2,Lc),Lc=1,3)
          call wegdmbs(xbs,ybs,bm,param,stk)
          bptilt=param(1)
	  bproll=param(2)
	  bptwst=param(3)
	  bpshft=param(5)	! Shear and Stagger swaped, Sep 18, 2004
	  bpslid=param(4)
	  bprise=param(6)
	return
	end
	subroutine errstop()
	write(6,1) 
c	write(6,*) 'NUPARM: Improper Usage'
c	write(6,*) 'Usage: NUPARM [-options] pdb_file_name'
	write(6,*) 'Options: -notdouble [default=double] '
	write(6,*) '     -notpdb [default=pdb]'
	write(6,*) '     -bpinf file_name_for_BasePairing_Information'
	write(6,*) '     -lsfit [default=cross-product method]'
	write(6,*) '     -c1c1 [default=C6-C8 line as Y-axis]'
	write(6,*) '     -cg [default=C6-C8 midpoint as BP center]'
	write(6,*) '     -single [default=Single Strand Parameter not',
     1 ' required'
	write(6,*) '     -orient "atom selection" [default=Not required]'
	write(6,*) '     -pp [default=Not required]'
	write(6,*) '     -torsion/-tor  [default=Not required]'
	write(6,*) '     -overlap/-ovl [default=Not required]'
	write(6,*) '     -axis [default=Not required]'
      write(6,*) '     -ww (to consider all base pairs as Watson-Crick 
     1type) [default=No]'
 1	format('NUPARM can be used in two ways:'/,
     1'  1. For a signle pdb file and'/,
     2'  2. For a trajectory file in PDB format as created by VMD (for e
     2xample).'/' In the first case:'/,'    NUPARM [-options] ',
     3 'pdb_file_name'/,' In the ssecond case:'/,
     4 '    NUPARM [-options] trj_file_name -trj -param name -output ',
     5 'out_file_name'/,' where the [-options] are as used in case 1.')
	stop
	return
	end

C***********************************************************************

!=======================================================================
!  CORRECTED SUBROUTINES : OVERLAP AND SURFACE
!
!=======================================================================
C***********************************************************************
      subroutine overlap(nb,ibapr,overlp)

      parameter ( natm = 999999 )
      parameter ( nrs = 9999 )
      include 'coordinates.h'
      include 'parameters.h'

      dimension ibapr(nrs),overlp(nrs)
      common/radii/surf(50000,3),rad(50),istp(50),ienp(50),nrad

      idpasq = 30     ! dots per angstrom square
c      npt=0
c      nrad=0
c      radl=0.0
c      open(unit=18,file='/usr/local/bin/surface.xyz',
c     1  STATUS='OLD')
c      do i=1,50000
c        npt=npt+1
c        read(18,25,end=311)radi,(surf(npt,l),l=1,3)
c        if(radi.ne.radl)then
c          nrad=nrad+1
c          rad(nrad)=radi
c          istp(nrad)=npt
c          radl=radi
c        endif
c        ienp(nrad)=npt
c      enddo
c311   continue
25    format(4f10.6)
C        write(99,*) 'nrad',nrad
C        write(99,25) rad(nrad),surf(nrad,1),surf(nrad,2),surf(nrad,3)
      ll=0
!======================= CHECK PRINT =============================
!      print*,nb,(nb-2)     ! Number of Bases
!=================================================================

      do inb=1,(nb-2),2
!      print*,'into the loop'
c      write(*,*)'done'

        nbp1 = ibapr(inb)
        nbp2 = ibapr(inb+1)
        nbp3 = ibapr(inb+2)
        nbp4 = ibapr(inb+3)

        call surface(nbp1,nbp2,0,0,0,nbpt1)
        call surface(nbp3,nbp4,0,0,0,nbpt2)
        call surface(nbp1,nbp2,nbp3,nbp4,1,nall)

        surf1=real(nbpt1)/real(idpasq)
        surf2=real(nbpt2)/real(idpasq)
        ovlp=real((nbpt1+nbpt2)-nall)
        ovlp=ovlp/(2.0*real(idpasq))
        ll = ll + 1
        overlp(ll)=ovlp
c        write(*,*)'Overlap = ',ovlp,' A^2'

!=================== CHECK PRINT ==============================
!      print*,nbpt1,nbpt2,nall,idpasq,ovlp
!==============================================================
      enddo
      close(unit=18)
      return 
      end

C***********************************************************************

C***********************************************************************
      subroutine surface(ibs1,ibs2,ibs3,ibs4,nst,ndots)

      parameter ( natm = 999999 )
      parameter ( nrs = 9999 )
      include 'coordinates.h'
      include 'parameters.h'

      dimension novl(100),rovl(100)
      common/radii/surf(50000,3),rad(50),istp(50),ienp(50),nrad
      common /trj/anstrj
      character*1 anstrj

      fac=0.2
      ir=0
      do n1=1,nseg
        if(((n1.eq.ibs1).or.(n1.eq.ibs2)).or.
     1  ((nst.eq.1).and.((n1.eq.ibs3).or.(n1.eq.ibs4))))then
          do n2=ists(n1),iens(n1)
            if((atnm(n2)(1:1).ne.'H').and.(atnm(n2)(3:3).ne.'''').and.
     1      (atnm(n2)(1:1).ne.'P').and.(atnm(n2)(3:3).ne.'P').and.
     2       (atnm(n2)(2:2).ne.'P'))then
              ir=ir+1
              novl(ir)=n2
              if(atnm(n2)(1:1).eq.'C')then
                atrad = 1.908
              elseif(atnm(n2)(1:1).eq.'N')then
                atrad = 1.824
              elseif(atnm(n2)(1:1).eq.'O')then
                atrad = 1.6612
c              elseif((atnm(n2).eq.' H21').or.(atnm(n2).eq.' H22').or.
c     1        (atnm(n2).eq.' H41').or.(atnm(n2).eq.' H42').or.
c     2        (atnm(n2).eq.' H61').or.(atnm(n2).eq.' H62').or.
c     3        (atnm(n2).eq.' H1 ').or.(atnm(n2).eq.' H3 '))then
c                atrad = 0.600
c              elseif((atnm(n2).eq.' H8 ').or.(atnm(n2).eq.' H2 '))then
c                atrad = 1.359
c              elseif(atnm(n2).eq.' H5 ')then
c                atrad = 1.459
c              elseif(atnm(n2).eq.' H6 ')then
c                atrad = 1.409
c              elseif((atnm(n2).eq.' H51').or.(atnm(n2).eq.' H52').or.
c     1        (atnm(n2).eq.' H53').or.(atnm(n2).eq.' H71').or.
c     2        (atnm(n2).eq.' H72').or.(atnm(n2).eq.' H73'))then
c                atrad = 1.487
              else
               atrad = 0.0
C 
C Added bhatta, Feb. 2015
C
              endif
              rovl(ir)=atrad + fac
            endif
          enddo
        endif
      enddo

      ndots=0
      do n3=1,ir
        x1=cr(novl(n3),1)
        y1=cr(novl(n3),2)
        z1=cr(novl(n3),3)
        r1=rovl(n3)
        do n4=1,nrad
          if(abs(r1-rad(n4)).le.0.0001)then
            do n5=istp(n4),ienp(n4)
              x11=x1+surf(n5,1)
              y11=y1+surf(n5,2)
              z11=z1+surf(n5,3)
              icnt=0
              do m=1,ir
                if(m.ne.n3)then
                  x2=cr(novl(m),1)
                  y2=cr(novl(m),2)
                  z2=cr(novl(m),3)
                  dist=sqrt((x11-x2)**2+(y11-y2)**2+(z11-z2)**2)
                  if(dist.le.rovl(m))then
                    icnt=1
                  endif
                endif
              enddo
              if(icnt.eq.0)then
                ndots=ndots+1
!=====================================================================
!      print surface (nst==0 for unstacked base pairs & nst==1 for 
!      stacked basepairs)
!=====================================================================
             if (nst ==1.and.anstrj.eq.'N')then
                write(37,591)'ATOM',ndots,atnm(novl(n3)),
     1          'NUC',ires(novl(n3)),x11,y11,z11,nst,rovl(n3)
             endif
!=====================================================================
              endif
            enddo
          endif
        enddo
      enddo
591   format(a4,2x,i5,1x,a4,1x,a3,3x,i3,4x,3f8.3,2x,i3,f8.2)
!==================== CHECK PRINT ============================================
!      print*,nst,icnt,ndots
!=============================================================================
      return
      end
C***********************************************************************
C***********************************************************************
      subroutine axisal(filep)
!this CODEX add basepair axis, base axis, and localhelix axis to the pdb
!files.  Author: Debasish Mukherjee, SINP 
	parameter ( nrs = 9999 )
      integer nres,natm,ixyzbps(30,nrs),countn
      real x0,x1,x2,x3,y0,y1,y2,y3,z0,z1,z2,z3
      character*40 filep
      character*200 fileo,xis,anme,file(3)
      character*80 deba
      character*5 debapart
      real xyzbse(4,30,30)

      iloop=2
C      write(*,*) 'the pdb file is',filep
!      goto 756
      open(unit=2,file=filep,status="old")                              !PDBfile
      nnf=index(filep,'.')
      fileo=filep
      fileo(nnf:nnf+11)='_axis_b.pdb'                                   !outfile name
      open(unit=3,file=fileo)                                           !outputfile
      fileo(nnf:nnf+12)='_axis_bp.pdb'
      open(unit=4,file=fileo)
      fileo(nnf:nnf+12)='_axis_lh.pdb'
      open(unit=7,file=fileo)
!-----------------------------------------------------------------------
      countn=0
      do while(1 .eq. 1)
      read(2,99,end=10) deba
99    format(a80)
      if(deba(1:4) .eq. "ATOM") then
      debapart(1:5)=deba(7:11)
c      read(debapart,'(i5)') natm                   !character to integer done
c      if(natm.eq.0) then
         countn=countn+1
         write(debapart,'(i5)') countn
         deba(7:11)=debapart(1:5)
c      endif
      read(deba(25:27),*) nres                      !character to integer done
      write(3,99) deba
      write(4,99) deba
      write(7,99) deba
      endif
      end do
10    continue
      natm=countn
      write(3,"(a3)")'TER'
      write(3,"(a3)")'END'
      write(4,"(a3)")'TER'
      write(4,"(a3)")'END'
      write(7,"(a3)")'TER'
      write(7,"(a3)")'END'
!-----------------------------------------------------------------------
      icnt=0
      icount=0
      rewind(unit=21)
      do while(1 .eq. 1)
      read(21,88,end=100) xis
88    format(a125)
!====== START OF ADD BASE AXIS========
      if(xis(1:7) .eq. "Base No")then
       icount=icount+1
       read(xis(10:13),*) icore
       if(icount .eq. 1)then
       item=icore
       else
       icount=0
       endif
       read(xis(21:30),*) x1
       read(xis(31:40),*) y1
       read(xis(41:50),*) z1
       read(xis(58:67),*) x2
       read(xis(68:77),*) y2
       read(xis(78:87),*) z2
       read(xis(95:104),*) x3
       read(xis(105:114),*) y3
       read(xis(115:124),*) z3
       xyzbse(1,icore,4)=x1
       xyzbse(1,icore,5)=y1
       xyzbse(1,icore,6)=z1
       xyzbse(1,icore,7)=x2
       xyzbse(1,icore,8)=y2
       xyzbse(1,icore,9)=z2
       xyzbse(1,icore,10)=x3
       xyzbse(1,icore,11)=y3
       xyzbse(1,icore,12)=z3
      endif
      if(xis(1:11) .eq. "Base ORIGIN")then
       icnt=icnt+1
       read(xis(12:19),*) x0
       read(xis(20:27),*) y0
       read(xis(28:35),*) z0
       if(icnt .eq. 1)then
       xyzbse(1,item,1)=x0
       xyzbse(1,item,2)=y0
       xyzbse(1,item,3)=z0
       else
       xyzbse(1,icore,1)=x0
       xyzbse(1,icore,2)=y0
       xyzbse(1,icore,3)=z0
       icnt=0
       endif
       endif
!====== END od ADD BASE AXIS=======
!====== START OF ADD BASE PAIR AXIS========
      if(xis(1:18) .eq. "BASE-PAIR ORIGIN-1")then
       read(xis(20:28),*) x0
       read(xis(30:38),*) y0
       read(xis(40:48),*) z0
       xyzbse(2,item-1,1)=x0
       xyzbse(2,item-1,2)=y0
       xyzbse(2,item-1,3)=z0
      endif
      if(xis(1:18) .eq. "BASE-PAIR ORIGIN-2")then
       read(xis(20:28),*) x0
       read(xis(30:38),*) y0
       read(xis(40:48),*) z0
       xyzbse(2,icore-1,1)=x0
       xyzbse(2,icore-1,2)=y0
       xyzbse(2,icore-1,3)=z0
      endif
      if(xis(1:11) .eq. "BP-XAXIS(1)")then
       read(xis(12:21),*) x1
       read(xis(22:31),*) y1
       read(xis(32:41),*) z1
       xyzbse(2,item-1,4)=x1
       xyzbse(2,item-1,5)=y1
       xyzbse(2,item-1,6)=z1
      endif
      if(xis(1:11) .eq. "BP-XAXIS(2)")then
       read(xis(12:21),*) x1
       read(xis(22:31),*) y1
       read(xis(32:41),*) z1
       xyzbse(2,icore-1,4)=x1
       xyzbse(2,icore-1,5)=y1
       xyzbse(2,icore-1,6)=z1
      endif
      if(xis(1:11) .eq. "BP-YAXIS(1)")then
       read(xis(12:21),*) x2
       read(xis(22:31),*) y2
       read(xis(32:41),*) z2
       xyzbse(2,item-1,7)=x2
       xyzbse(2,item-1,8)=y2
       xyzbse(2,item-1,9)=z2
      endif
      if(xis(1:11) .eq. "BP-YAXIS(2)")then
       read(xis(12:21),*) x2
       read(xis(22:31),*) y2
       read(xis(32:41),*) z2
       xyzbse(2,icore-1,7)=x2
       xyzbse(2,icore-1,8)=y2
       xyzbse(2,icore-1,9)=z2
      endif
      if(xis(1:11) .eq. "BP-ZAXIS(1)")then
       read(xis(12:21),*) x3
       read(xis(22:31),*) y3
       read(xis(32:41),*) z3
       xyzbse(2,item-1,10)=x3
       xyzbse(2,item-1,11)=y3
       xyzbse(2,item-1,12)=z3
      endif
      if(xis(1:11) .eq. "BP-ZAXIS(2)")then
       read(xis(12:21),*) x3
       read(xis(22:31),*) y3
       read(xis(32:41),*) z3
       xyzbse(2,icore-1,10)=x3
       xyzbse(2,icore-1,11)=y3
       xyzbse(2,icore-1,12)=z3
      endif
!====== END of ADD BASE PAIR AXIS=======
!====== START OF LOCAL AXIS=============
      if(xis(12:27) .eq. "Local Helix Axis")then
       if(xis(35:37) .ne."nan")then
       read(xis(28:37),*) x3
          read(xis(38:47),*) y3
          read(xis(48:57),*) z3
          xyzbse(3,item-1,10)=x3
          xyzbse(3,item-1,11)=y3
          xyzbse(3,item-1,12)=z3
       endif
      endif
      if(xis(1:18) .eq. "Local Helix Origin")then
      iloop=iloop+2
       if(xis(26:28) .ne."nan")then
          read(xis(20:28),*) x0
          read(xis(30:38),*) y0
          read(xis(40:48),*) z0
          xyzbse(3,item-1,1)=x0
          xyzbse(3,item-1,2)=y0
          xyzbse(3,item-1,3)=z0
       endif
      endif
!====== END OF LOCAL AXIS=============
      end do
100   continue
!-----------------------------------------------------------------------
      ilooph=iloop/2
      DO l=1,iloop
       DO m=1,3
       xyzbse(1,l,m)=xyzbse(1,l,m)
       xyzbse(1,l,m+3)=xyzbse(1,l,m)+xyzbse(1,l,m+3)*2
       xyzbse(1,l,m+6)=xyzbse(1,l,m)+xyzbse(1,l,m+6)*2
       xyzbse(1,l,m+9)=xyzbse(1,l,m)+xyzbse(1,l,m+9)*2
       xyzbse(2,l,m)=xyzbse(2,l,m)
       xyzbse(2,l,m+3)=xyzbse(2,l,m)+xyzbse(2,l,m+3)*4
       xyzbse(2,l,m+6)=xyzbse(2,l,m)+xyzbse(2,l,m+6)*8
       xyzbse(2,l,m+9)=xyzbse(2,l,m)+xyzbse(2,l,m+9)*2
       xyzbse(3,l,m)=xyzbse(3,l,m)
       xyzbse(3,l,m+9)=xyzbse(3,l,m)+xyzbse(3,l,m+9)*2
!       write(*,*) n,l,m
       end do
      end do
      ip=2
!      Do j=1,1
       DO n=1,10,3
         nres=nres+1
         DO L=1,iloop
           natm=natm+1
           if(n .eq. 1)then
           anme='O'
           endif
           if(n .eq. 4)then
             anme='X'
           endif
           if(n .eq. 7)then
           anme='Y'
           endif
           if(n .eq. 10)then
           anme='Z'
           endif
           ixyzbps(n,L)=natm
           write(3,94) 'HETATM',natm,anme,'BPR','A',nres,(xyzbse(1,L,I),
     1 I=n,n+2)
           if(L .le. ilooph)then
           write(4,94) 'HETATM',natm,anme,'BPR','A',nres,(xyzbse(2,L,I),
     1 I=n,n+2)
           endif
           if(L .lt. ilooph)then
             if((n .lt. 4) .or. (n .gt. 9))then
             write(7,94) 'HETATM',natm,anme,'BPR','A',nres,(xyzbse(3,L,
     1 I),I=n,n+2)
             endif
           endif
94           format (a6,2x,i3,2x,a1,3x,a3,1x,1a,2x,i2,4x,3F8.3)
         end do
       end do
      write(3,"(a3)")'TER'
      write(4,"(a3)")'TER'
      write(7,"(a3)")'TER'
      DO ij=4,10,3
       DO io=1,iloop
       write(3,"(a6,2x,i3,2x,i3)") 'CONECT',ixyzbps(1,io),ixyzbps(ij,io)
       if(io .lt. ilooph)then
       if((ij .lt. 4) .or. (ij .gt. 9))then
       write(7,"(a6,2x,i3,2x,i3)") 'CONECT',ixyzbps(1,io),ixyzbps(ij,io)
         endif
       endif
       ENDDO
      ENDDO
      DO ij=4,10,3
       DO io=1,ilooph
       write(4,"(a6,2x,i3,2x,i3)") 'CONECT',ixyzbps(1,io),ixyzbps(ij,io)
       ENDDO
      ENDDO
      write(3,"(a3)")'END'
      write(4,"(a3)")'END'
      write(7,"(a3)")'END'
      write(*,*) 'Number of residue',iloop
      close(unit=21)
      close(unit=2)
      close(unit=3)
      close(unit=4)
      close(unit=7)
!756   continue
      return
      end
      function distn(xx,yy)
      dimension xx(3),yy(3)
      DIST=0.0
      DO 2006 K=1,3
         DIST=DIST+(Xx(K)-Yy(K))*(Xx(K)-Yy(K))
 2006 CONTINUE
      distn=sqrtf(dist)
      return
      end

C----------------------------------------------------------------------
C  SUBROUTINE TO FIND OUT DIRECTION TOWARD MAJOR GROOVE, TO FIX X-AXIS.
C----------------------------------------------------------------------
      SUBROUTINE DIRECX(AXIS,I,possblz,foundz1,foundz2,foundz3)
      DOUBLE PRECISION AXIS,possblz,foundz1,foundz2,foundz3
      DIMENSION AXIS(2,3),possblz(3)

      IND = 0
      ANG =  possblz(1)*foundz1
      ang = ang + possblz(2)*foundz2
      ang = ang + possblz(3)*foundz3
      
      IF(ANG.LT.0.0) THEN
           DO 2002 KJ=1,3
                      AXIS(I,KJ) = -AXIS(I,KJ)
 2002      CONTINUE
      ENDIF
      RETURN
      END
      subroutine readsurf
      common/radii/surf(50000,3),rad(50),istp(50),ienp(50),nrad

      write(99,*) 'Overlap calculation requested'
      open(unit=18,file='/usr/local/bin/surface.xyz',
     1  STATUS='OLD')
      npt=0
      nrad=0
      radl=0.0
      do i=1,50000
        npt=npt+1
        read(18,25,end=311)radi,(surf(npt,l),l=1,3)
        if(radi.ne.radl)then
          nrad=nrad+1
          rad(nrad)=radi
          istp(nrad)=npt
          radl=radi
        endif
        ienp(nrad)=npt
      enddo
311   continue
25    format(4f10.6)
        return
        end
