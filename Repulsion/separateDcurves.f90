    MODULE SDC_globals
    integer resFileNo, freqsNum, SDcurvesNum
    character (len=40) curvesDir
    real*8 fOld, alfaOld, SDCstep, SDCeps, fMax
    real*8 fMin, fStep, dzMin, dzMax, haminStep, haminEps
    complex*16 pOld, p_c
    
    real*8, allocatable:: freqs(:), dzetas(:), resK11(:), p1(:,:), p2(:,:)
    namelist /basicInfo/SDcurvesNum, fmax, SDCstep, SDCeps, curvesDir
    namelist /startPoints/ p1, p2
    namelist /simpleDC/ freqsNum, fMin, dzMin, dzMax, haminStep, haminEps
    
    CONTAINS 
        subroutine initSeparateDcurves
        IMPLICIT NONE
            open(unit=1,file='SeparateDcurves settings.txt',status='old')
            read(1, basicInfo)
            allocate(p1(3,SDcurvesNum), p2(3,SDcurvesNum))
            read(1, startPoints)
            read(1, simpleDC)
            close(1)       
        end
    END MODULE
    
    
    
    SUBROUTINE RPolesDots
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
        write(301,'(A)') "%  f, MHz;           Re(dzeta);"
        fstep = (fmax - fmin)/freqsNum;
        print*, 'Regular, non separate d curver plot has been started'
        do i = 1, freqsNum
            f = fmin + fstep*(i-1); w = f*2d0*pi; kap2=w*hs(1)/vs;
            call Hamin(RDabs, dzMin, dzMax, haminStep, haminEps, 20, dz, Ndz)
            do j = 1, Ndz
                write(301, 1) f, dz(j)*w, dz(j)
            enddo    
        enddo
        close(301)
        1       format(3F30.16)
    END SUBROUTINE RPolesDots
    
    
                                                                ! auxiliary function for real poles finding by Hamin
    real(8) function arcRDabs(angle)
    use Mult_Glob; 
    use SDC_globals;
    implicit none
	real(8) gm,angle; 
    complex*16 alf,det
	common/gm/gm/det/det
        w = 2d0*pi*(cos(angle)*SDCstep + fOld); f=cos(angle)*SDCstep + fOld;
	    alf = (sin(angle)*SDCstep + alfaOld)*(1d0,0d0)
	    call MultiK_An(alf,gm,0d0)
        arcRDabs = 1d0/(abs(Kaz(1,1))+abs(Kaz(2,2))+abs(Kaz(3,3)))
    end    
    
    
    function arcHaminDelta(angle) result(delta)
    use Mult_Glob; 
    use SDC_globals;
    implicit none
	real*8 gm,angle, delta
    complex*16 alf,det, p_cur
	common/gm/gm/det/det
        p_cur = pOld + p_c*exp(cci*angle)
        f = real(p_cur); w = 2d0*pi*f;
        alf = imag(p_cur)*(1d0,0d0);
	    call MultiK_An(alf,gm,0d0)
        delta = 1d0/(abs(Kaz(1,1))+abs(Kaz(2,2))+abs(Kaz(3,3)))
    end
    
    


    
    
    subroutine RPoleTracer(p1, p2, fileNo, arcStep, radStep, fMax, logNo)
    use SDC_globals, only: pOld, p_c, SDCeps
    use Mult_Glob
    implicit none
    logical halfed
    integer fileNo, logNo, Ndz, forward
    real*8 arcStep, radStep, fMax, dz(20), arcHaminDelta, oldRadStep
    real*8, allocatable:: abs_dz(:)
    complex*16 p, p1, p2
    external arcHaminDelta
        write(fileNo, 2) p1
        write(fileNo, 2) p2
        oldRadStep = RadStep
        do 
            p_c = p2 - p1; p_c = p_c/abs(p_c)*radStep;
            pOld = p2;
            if (real(p2) > fmax .OR. real(p2)<=0d0) exit;
            call Hamin(arcHaminDelta, -0.75d0*pi, 0.75d0*pi, arcStep, SDCeps, 20, dz, Ndz)
            if (Ndz == 1) then 
                p1 = p2; p2 = p2 + p_c*exp(cci*dz(1));
                write(fileNo, 2) real(p2), imag(p2), 0d0
                radStep = oldRadStep;
            else if (Ndz>1) then
                !print*, 'RPolesTracer: To many poles!', Ndz
                allocate(abs_dz(Ndz))
                abs_dz = abs(dz(1:Ndz))
                forward = minloc(abs_dz, DIM = 1)
                deallocate(abs_dz)
                p1 = p2; p2 = p2 + p_c*exp(cci*dz(forward));
                write(fileNo, 2) real(p2), imag(p2), 0d0
                radStep = oldRadStep;
            else 
                !print*, 'RPolesTracer: No poles at all!', Ndz
                radStep = radStep/1.5d0; halfed = .TRUE.;
                !print*, 'RPolesTracer: Half step!', radStep
            endif
        enddo
 2       format(3F30.16)   
    end
    
    
    subroutine RPoleCurves ! открывает папку с файлами дисп. кривых, строит кривые от начальных точек f_sp(i), alfa_sp(i)
    use SDC_globals;
    implicit none
    integer i, fileNum
    complex*16 theP1, theP2
    character(len=20) fileName, str
    external str
        call initSeparateDcurves
        print*, 'Separate disp curves plotting has been started!'
        do i = 3, SDcurvesNum
            fileName = str(i)
            fileNum = i+300
            theP1 = p1(1, i) + (0d0, 1d0)*p1(2, i); 
            theP2 = p2(1, i) + (0d0, 1d0)*p2(2, i); 
            print*, 'Curve ', i, ' plotting has been started!'
            open(unit=fileNum, file='.\DataFigs\separateDcurves\'//trim(curvesDir)//'\Dcurves\'//trim(fileName)//'.txt', FORM='FORMATTED')
            write(fileNum,'(A)') "%  f, MHz;           Re(dzeta);              Im(dzeta)"
            call RPoleTracer(theP1, theP2, fileNum, 1d-2, SDCstep, fmax, 2)
            close(fileNum)
            print*, 'Curve ', i, 'is done!'
        enddo    
    end subroutine RPoleCurves
    
    
    
    
    
    
    
    subroutine resFromRPoles ! открывает папку с файлами дисп. кривых, строит соотв. файлы кривых амплитуд
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
            write(newFileNum,'(A)') "%  f, MHz;           Re(dzeta);              abs(K11), abs(K13), abs(K31), abs(K33)"
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
    
    !
    !subroutine SCP1(f0, alfa0, file)  ! Ќаходит точки дисперсионной кривой и записывает в файл соотв значени€ частоты и полюса, не переписывалась еще с момента как заработала
    !use SDC_globals;
    !use Mult_Glob; 
    !implicit none
    !integer  Ndz, choice, iterno, j, file
    !real*8 f0, alfa0, fNew, alfaNew, step, dz(20), psi, arcRDabs
    !external arcRDabs
    !!    !                                                                    первые шаги, подготовка к автоматике
    !    step = SDCstep; fOld = f0; alfaOld = alfa0;
    !    call Hamin(arcRDabs, 0d0, pi, 2d-3, SDCeps, 20, dz, Ndz)
    !    alfaNew = sin(dz(1))*SDCstep + alfaOld; fNew = cos(dz(1))*SDCstep + fOld;
    !    psi = atan( (fNew-fOld)/(alfaNew-alfaOld) );      
    !    w = 2d0*pi*fNew       
    !    write(file, '(4E15.6E3)') fNew, alfaNew
    !    fNew = fOld; alfaNew = alfaOld;
    !!                                                                         автоматический режим
    !    do 
    !        iterno= 0;
    !        do
    !            call Hamin(arcRDabs, pi/3d0-psi, 2d0*pi/3d0-psi, 0.5d-3, SDCeps, 4, dz, Ndz)
    !            if (Ndz>1) then
    !                !print*, Ndz
    !                choice = 1; 
    !                do j = 2, Ndz
    !                    if ( abs(dz(j)-(pi/2d0 - psi)) < abs(dz(j-1)-(pi/2d0 - psi)) ) choice = j
    !                enddo
    !                exit;
    !            else
    !                if (step < SDCstep) SDCstep = SDCstep*2d0
    !                 choice = 1; exit;
    !            endif    
    !            iterno = iterno + 1; if (iterno > 5) exit;
    !        enddo
    !        alfaNew = sin(dz(choice))*SDCstep + alfaOld; fNew = cos(dz(choice))*SDCstep + fOld;
    !        w = 2d0*pi*fNew        
    !        write(file, '(4E15.6E3)') fNew, alfaNew
    !        psi = atan( (fNew-fOld)/(alfaNew-alfaOld) ); fOld = fNew; alfaOld = alfaNew;
    !        if (fNew>fmax) exit;
    !    enddo 
    !end subroutine SCP1
    !
    
    !subroutine plotAllcurves ! открывает папку с файлами дисп. кривых, строит кривые от начальных точек f_sp(i), alfa_sp(i)
    !use SDC_globals;
    !implicit none
    !integer i, fileNum
    !complex*16 p1, p2
    !character(len=20) fileName, str
    !external str
    !    call initSeparateDcurves
    !    print*, 'Separate disp curves plotting has been started!'
    !    do i = 1, SDcurvesNum
    !        fileName = str(i)
    !        fileNum = i+300
    !        print*, 'Curve ', i, ' plotting has been started!'
    !        open(unit=fileNum, file='.\DataFigs\separateDcurves\'//trim(curvesDir)//'\Dcurves\'//trim(fileName)//'.txt', FORM='FORMATTED')
    !        call SCP1(f_sp(i), alfa_sp(i), fileNum)
    !        close(fileNum)
    !        print*, 'Curve ', i, 'is done!'
    !    enddo    
    !end subroutine plotAllcurves