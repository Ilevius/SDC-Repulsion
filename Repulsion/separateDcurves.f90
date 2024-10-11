    MODULE SDC_globals
    integer resFileNo, freqsNum, SDcurvesNum, SDcurvesNumComplex, dotsNum, mx, IPP
    character (len=40) curvesDir
    real*8 fOld, alfaOld, SDCstep, SDCstepComplex, SDCeps, SDCepsComplex, fMax, fmaxComplex
    real*8 fMin, fStep, dzMin, dzMax, haminStep, haminEps
    complex*16 dzStep
    
    
    real*8, allocatable:: freqs(:), dzetas(:), resK11(:), f_sp(:), alfa_sp(:), p1(:,:), p2(:,:)
    namelist /basicInfo/SDcurvesNum, fmax, SDCstep, SDCeps, curvesDir
    namelist /startPoints/ f_sp, alfa_sp
    namelist /simpleDC/ dotsNum, fMin, dzMin, dzMax, haminStep, haminEps
    namelist /basicInfoComplex/SDcurvesNumComplex, fmaxComplex, SDCstepComplex, dzStep, SDCepsComplex, mx, IPP
    namelist /startPointsComplex/ p1, p2
    
    
    CONTAINS 
        subroutine initSeparateDcurves
        IMPLICIT NONE
            open(unit=1,file='SeparateDcurves settings.txt',status='old')
            read(1, basicInfo)
            allocate(f_sp(SDcurvesNum), alfa_sp(SDcurvesNum))
            read(1, startPoints)
            read(1, simpleDC)
            read(1, basicInfoComplex)
            allocate(p1(3, SDcurvesNumComplex), p2(3, SDcurvesNumComplex))
            read(1, startPointsComplex)
            close(1)       
        end
    END MODULE
    
    
    
    SUBROUTINE simpleDcurves
    use Mult_Glob;            ! GIVES US w, f, pi
    use SDC_globals;
    IMPLICIT NONE
    real*8 dz(20)
    integer i, j, Ndz
    real*8 RDabs, RDabs_alf, vs, kap2
    external RDabs,RDabs_alf
        call initSeparateDcurves
        dz = 0d0;
        vs = sqrt(C_6(6,6,1)/rho(1))
        open(unit=301, file='.\DataFigs\separateDcurves\'//trim(curvesDir)//'\Dcurves\simpleDcurves.txt', status='unknown')
        fstep = (fmax - fmin)/dotsNum;
        print*, 'Regular, non separate d curver plot has been started'
        do i = 1, dotsNum
            f = fmin + fstep*(i-1); w = f*2d0*pi; kap2=w*hs(1)/vs;
            call Hamin(RDabs, dzMin, dzMax, haminStep, haminEps, 20, dz, Ndz)
            do j = 1, Ndz
                write(301, '(5E15.6E3)') f, dz(j)*w, dz(j)
            enddo    
        enddo
        close(301)
    END SUBROUTINE simpleDcurves
    
    
                                                                ! auxiliary function for real poles finding by Hamin
    real(8) function arcRDabs(angle)
    use Mult_Glob; 
    use SDC_globals;
    implicit none
	real(8) gm,angle; 
    complex*16 alf,c,det
	common/gm/gm/det/det
        w = 2d0*pi*(cos(angle)*SDCstep + fOld)
	    alf = (sin(angle)*SDCstep + alfaOld)*(1d0,0d0)
	    call MultiK_An(alf,gm,0d0)
        c = 1d0/Kaz(3,3)
	    arcRDabs = abs(c)
    end
    
    
    
    function CrootDeltaPLS(alf) result(c) 
    use Mult_Glob; 
    implicit none
    real(8) gm,angle; 
    complex*16 alf,c,det
	common/gm/gm/det/det
        call MultiK_An(alf,gm,0d0)
        c = 1d0/Kaz(3,3)
    end
    
    
    function testCrootDelta(alfa) result(res)
    use Mult_Glob
    implicit none
    complex*16 alfa, res
        res = alfa - sqrt(f+0.001) - cci*sqrt(f)
    end
    

    subroutine SCP1(f0, alfa0, file)  ! Ќаходит точки дисперсионной кривой и записывает в файл соотв значени€ частоты и полюса, не переписывалась еще с момента как заработала
    use SDC_globals;
    use Mult_Glob; 
    implicit none
    integer  Ndz, choice, iterno, j, file
    real*8 f0, alfa0, fNew, alfaNew, step, dz(20), psi, arcRDabs
    external arcRDabs
    !    !                                                                    первые шаги, подготовка к автоматике
        step = SDCstep; fOld = f0; alfaOld = alfa0;
        call Hamin(arcRDabs, 0d0, pi, 1d-3, SDCeps, 20, dz, Ndz)
        alfaNew = sin(dz(1))*SDCstep + alfaOld; fNew = cos(dz(1))*SDCstep + fOld;
        psi = atan( (fNew-fOld)/(alfaNew-alfaOld) );      
        w = 2d0*pi*fNew       
        write(file, '(4E15.6E3)') fNew, alfaNew
        fNew = fOld; alfaNew = alfaOld;
    !                                                                         автоматический режим
        do 
            iterno= 0;
            do
                call Hamin(arcRDabs, pi/3d0-psi, 2d0*pi/3d0-psi, 0.5d-3, SDCeps, 4, dz, Ndz)
                if (Ndz>1) then
                    !print*, Ndz
                    choice = 1; 
                    do j = 2, Ndz
                        if ( abs(dz(j)-(pi/2d0 - psi)) < abs(dz(j-1)-(pi/2d0 - psi)) ) choice = j
                    enddo
                    exit;
                else
                    if (step < SDCstep) SDCstep = SDCstep*2d0
                     choice = 1; exit;
                endif    
                iterno = iterno + 1; if (iterno > 5) exit;
            enddo
            alfaNew = sin(dz(choice))*SDCstep + alfaOld; fNew = cos(dz(choice))*SDCstep + fOld;
            w = 2d0*pi*fNew        
            write(file, '(4E15.6E3)') fNew, alfaNew
            psi = atan( (fNew-fOld)/(alfaNew-alfaOld) ); fOld = fNew; alfaOld = alfaNew;
            if (fNew>fmax) exit;
        enddo 
    end subroutine SCP1
    
    
    
    subroutine plotAllcurves ! открывает папку с файлами дисп. кривых, строит кривые от начальных точек f_sp(i), alfa_sp(i)
    use SDC_globals;
    implicit none
    integer i, fileNum
    character(len=20) fileName, str
    external str
        call initSeparateDcurves
        print*, 'Separate disp curves plotting has been started!'
        do i = 1, SDcurvesNum
            fileName = str(i)
            fileNum = i+300
            print*, 'Curve ', i, ' plotting has been started!'
            open(unit=fileNum, file='.\DataFigs\separateDcurves\'//trim(curvesDir)//'\Dcurves\'//trim(fileName)//'.txt', FORM='FORMATTED')
            call SCP1(f_sp(i), alfa_sp(i), fileNum)
            close(fileNum)
            print*, 'Curve ', i, 'is done!'
        enddo    
    end subroutine plotAllcurves
    
    
    
    
    
    
    
    subroutine residuesFromDCurves ! открывает папку с файлами дисп. кривых, строит соотв. файлы кривых амплитуд
    use SDC_globals;
    implicit none
    integer i, j, oldFileNum, newFileNum, io
    character(len=20) fileName, str
    real*8 freq, dzeta
    complex*16 K11, K13, K31, K33
    external str
        call initSeparateDcurves
    ! главный цикл пробегает по файлам дисп кривых, открывает их и читает
        do i = 1, SDcurvesNum
            fileName = str(i); oldfileNum = i + 300; newFileNum = i + 400;
            print*, 'Curve ', i, 'residues plotting has been started!'
            open(unit=oldfileNum, file='.\DataFigs\separateDcurves\'//trim(curvesDir)//'\Dcurves\'//trim(fileName)//'.txt', status='old')
            open(unit=newFileNum, file='.\DataFigs\separateDcurves\'//trim(curvesDir)//'\Residues\'//trim(fileName)//'.txt', FORM='FORMATTED')
            do
                ! цикл читает файл дисп кривой построчно, находит вычеты и записывает в строку соответствующего файла вычетов
                read(oldFileNum, *, iostat=io) freq, dzeta
                if (io/=0) exit
                call residueAtAdzeta(freq, dzeta, K11, K13, K31, K33)
                write(newFileNum, '(7E15.6E3)') freq, dzeta, abs(K11), abs(K13), abs(K31), abs(K33)
            enddo
            close(newFileNum); close(oldFileNum);
            print*, 'Curve ', i, 'is done!'
        enddo
    end     
    
    
    subroutine residueAtAdzeta(freq, dzeta, K11, K13, K31, K33) ! ¬ычисл€ет некоторые вычеты в частоте freq и полюсе dzeta
    use SDC_globals
    use Mult_Glob
    implicit none
    real*8 freq, dzeta, hres, gm
    complex*16 K11, K13, K31, K33
    complex(8) alf
        w = 2d0*pi*freq
        hres = 1d3*SDCeps
        if((dzeta-hres) < 0d0) dzeta = dzeta-eps
        
        alf = dzeta + hres	
        gm = 0d0
        call MultiK_An(alf,gm,0d0); 
        
        K11 = Kaz(1,1); K31 = Kaz(3,1)  
        K13 = Kaz(1,3); K33 = Kaz(3,3)  
           
        alf = dzeta - hres	
        call MultiK_An(alf,gm,0d0); 
        
        K11 = hres*(K11-Kaz(1,1))/2d0; 
        K31 = hres*(K31-Kaz(3,1))/2d0  
        K13 = hres*(K13-Kaz(1,3))/2d0 
        K33 = hres*(K33-Kaz(3,3))/2d0  
    end
    
    
    subroutine countFileLines(fileNo, nlines)
    implicit none
    integer fileNo, nlines, io
        nlines = 0
        DO
          READ(fileNo,*,iostat=io)
          IF (io/=0) EXIT
          nlines = nlines + 1
        END DO
    end
    
    
    character(len=20) function str(k)
!   "Convert an integer to string."
    integer, intent(in) :: k
        write (str, *) k
        str = adjustl(str)
    end function str
    
    
   
!                                                         _______  _______  _______  _______ _________
!                                                        (  ____ \(  ____ )(  ___  )(  ___  )\__   __/
!                                                        | (    \/| (    )|| (   ) || (   ) |   ) (   
!                                                        | |      | (____)|| |   | || |   | |   | |   
!                                                        | |      |     __)| |   | || |   | |   | |   
!                                                        | |      | (\ (   | |   | || |   | |   | |   
!                                                        | (____/\| ) \ \__| (___) || (___) |   | |   
!                                                        (_______/|/   \__/(_______)(_______)   )_(   
                                             
    
    subroutine Croot13(oldP)
    use Mult_Glob
    use SDC_globals
    implicit none
    real*8 oldP(3), newP(3)
    complex*16 theRoot, step, CrootDeltaPLS, testCrootDelta
    external CrootDeltaPLS, testCrootDelta
        w = 2d0*pi*oldP(3); f = oldP(3);
        theRoot = oldP(1)+(0d0, 1d0)*oldP(2)
        call CROOTW25(theRoot,dzStep,SDCepsComplex,mx,IPP,CrootDeltaPLS)
        newP(1) = real(theRoot); newP(2) = imag(theRoot); newP(3) = oldP(3);
        oldP = newP
    end
    
    
    subroutine CRcurve(sp1, sp2, fileNum)
    use SDC_globals
    implicit none
    integer fileNum
    real*8 sp1(3), sp2(3), sp3(3), n(3)
        write(fileNum, '(5E15.6E3)') p1
        write(fileNum, '(5E15.6E3)') p2
        do 
            n = sp2-sp1
            n = n/sqrt(n(1)**2+n(2)**2+n(3)**2)
            sp3 = sp2 + n*SDCstepComplex
            call Croot13(sp3)
            write(fileNum, '(5E15.6E3)') sp3
            sp1 = sp2; sp2 = sp3;
            if(sp3(3)>fmaxComplex .or. sp3(3)<0d0) exit
        enddo
    end
        
    

    subroutine PlotAllCRcurves ! открывает папку с файлами дисп. кривых, строит кривые от начальных точек f_sp(i), alfa_sp(i)
    use SDC_globals;
    implicit none
    integer i, fileNum
    character(len=20) fileName, str
    external str
        call initSeparateDcurves
        print*, 'Separate complex disp curves plotting has been started!'
        do i = 1, SDcurvesNum
            fileName = str(i)
            fileNum = i+300
            print*, 'Curve ', i, ' plotting has been started!'
            open(unit=fileNum, file='.\DataFigs\separateDcurves\'//trim(curvesDir)//'\CDcurves\'//trim(fileName)//'.txt', FORM='FORMATTED')
            call CRcurve(p1(:,i), p2(:,i), fileNum)
            close(fileNum)
            print*, 'Curve ', i, 'is done!'
        enddo    
    end subroutine PlotAllCRcurves