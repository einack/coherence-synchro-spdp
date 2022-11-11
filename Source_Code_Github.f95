PROGRAM Network
    IMPLICIT NONE
    
    ! ******************************************************************************
    ! Parameters for network topology
    ! ******************************************************************************
    INTEGER :: ind, iiv, inode, nodes, neigh 
    DOUBLE PRECISION :: rnd 
    INTEGER, ALLOCATABLE :: unwired_connections(:), weights(:,:), c1(:)
    ! ******************************************************************************
    DOUBLE PRECISION, ALLOCATABLE :: S1(:,:), V(:,:), m(:,:), h(:,:), n(:,:), u(:,:)
    DOUBLE PRECISION :: pi, t, ntk 
    INTEGER, PARAMETER::  trans_time = 75000  !Transient time 
    
    ! ******************************************************************************
      ! Variables for calculating the coefficient of variation 
      ! ******************************************************************************
    DOUBLE PRECISION, ALLOCATABLE :: S(:,:), ISI1(:,:) 
    DOUBLE PRECISION, ALLOCATABLE :: SUM_ISI1(:), SUM_ISI_sq1(:)
    DOUBLE PRECISION, ALLOCATABLE :: MEAN_ISI1(:), MEAN_ISI_sq1(:)
    DOUBLE PRECISION :: SUM_MEAN_ISI1, SUM_MEAN_ISI_sq1
    DOUBLE PRECISION :: MEAN1, MEAN_sq1
    INTEGER :: r
    ! ******************************************************************************
  
    DOUBLE PRECISION :: v1, v2, v3, v4
    DOUBLE PRECISION :: m1, m2, m3, m4
    DOUBLE PRECISION :: h1, h2, h3, h4
    DOUBLE PRECISION :: n1, n2, n3, n4
    DOUBLE PRECISION :: p1, p2, p3, p4
  
    DOUBLE PRECISION :: noise_V, noise_m, noise_h, noise_n
    DOUBLE PRECISION :: alpha_m, alpha_h, alpha_n
    DOUBLE PRECISION :: beta_m, beta_h, beta_n
    
    DOUBLE PRECISION :: a1, a2, a3, a4, a5, a6, a7, a8
    DOUBLE PRECISION :: a11, a22, a33, a44
    
    ! ******************************************************************************
    ! Parameters for the HH model
    ! ******************************************************************************
    DOUBLE PRECISION, PARAMETER:: g_Na = 120.0d0, g_K = 36.0d0, g_L = 0.3d0
    DOUBLE PRECISION, PARAMETER:: V_Na = 50.d0, V_K = -77.d0, V_L = -54.4d0
    DOUBLE PRECISION, PARAMETER:: C_m = 1.0d0 , I_0 = 6.0d0
    DOUBLE PRECISION, PARAMETER:: rho_Na = 60.0d0, rho_K = 18.0d0 
    ! ******************************************************************************
    
    ! ******************************************************************************
    ! Parameters for the synapses
    ! ******************************************************************************
    !DOUBLE PRECISION ::  gstar!, test_time
    
    DOUBLE PRECISION, ALLOCATABLE:: dS1(:) 
    DOUBLE PRECISION, ALLOCATABLE:: W1(:,:), dW1(:)
    DOUBLE PRECISION, PARAMETER:: min_W1 = 0.0001d0, max_W1 = 1.0d0
    DOUBLE PRECISION, PARAMETER:: mean_W1 = 0.1d0 
    DOUBLE PRECISION, PARAMETER:: sd_W1 = 0.02d0  
    
    DOUBLE PRECISION, PARAMETER:: learning_rate = 0.0001d0
    DOUBLE PRECISION, PARAMETER:: V_r = -75.0d0 !, tau_d = 5.0d0, tau_r = 0.5d0
    DOUBLE PRECISION, PARAMETER:: tau_1 = 20.0d0, tau_2 = 20.0d0
    DOUBLE PRECISION, PARAMETER:: A_1 = 1.0d0 , A_2 = 0.5d0 
    DOUBLE PRECISION, PARAMETER:: B_1 = 1.0d0 , B_2 = 1.1d0
    ! *****************************************************************************
    
    INTEGER, ALLOCATABLE :: SEED(:)
    INTEGER:: iv, k, p, j, c, Nstep, iseed, q  !, c2
    DOUBLE PRECISION, PARAMETER:: thres = 20.0d0, dt = 0.01d0
    INTEGER, PARAMETER:: realizations = 50 , noise_data = 34
    DOUBLE PRECISION, PARAMETER:: kappa1 = 0.25d0, tau = 5.0d0
    DOUBLE PRECISION, PARAMETER:: sigma = 0.0d0, tol = 0.01d0
    DOUBLE PRECISION :: start, finish, tic, tac !, pcum
    DOUBLE PRECISION :: CV(realizations), CV_sum, CV_sum2, CV_ave
    DOUBLE PRECISION :: CV_sd
    DOUBLE PRECISION :: Sum_1(realizations)
    DOUBLE PRECISION :: s_deviation(realizations)
    DOUBLE PRECISION :: Sum_2, s_deviation_2 
    DOUBLE PRECISION :: avg_w, avg_sd
    DOUBLE PRECISION :: beta, channel_noise(noise_data)
  
  ! ******************************************************************************
  ! Synchro index and order parameter
   ! ********************************************************************************
    DOUBLE PRECISION, ALLOCATABLE :: sum_network_V(:), sum_network2_V(:), std_dev(:)
    DOUBLE PRECISION :: mean_network_V, mean_network2_V
    DOUBLE PRECISION :: sum_std_dev, mean_std_dev(realizations) 
    DOUBLE PRECISION :: syn_index_sum, syn_index_sum2, synronization, synronization_sd
    ! ******************************************************************************
  
  
    call cpu_time(start)
    
    OPEN(11,File="Average_synaptic_weight.txt",Status='unknown')
    OPEN(12,File="Coefficient_of_variation.txt",Status='unknown')
    OPEN(13,File="Synchronization_index.txt",Status='unknown')
    OPEN(10,File="noise_data.txt", Status='unknown', action='read')
    !OPEN(14,File="small_world.txt",Status='unknown')


    pi = 4.0*ATAN(1.0)

    
     ! Network parameters initialization
      OPEN(UNIT=17,FILE='parameters_sw.txt', STATUS='unknown')
      READ(17,*)
      READ(17,*) nodes, neigh, beta, Nstep
      CLOSE(17)
      
      PRINT*, "Number of nodes:", nodes
      PRINT*, "Number of neighbours:", neigh
      PRINT*, "Rewiring probability:", beta
      PRINT*, "Integration time steps:", Nstep
      PRINT*
      
      IF(neigh > NINT(nodes/2.0)-1)THEN
        PRINT*, 'Error: The number of neighbours cannot exceed half of number of nodes'
        STOP 
      END IF
      
      IF(beta > 1.0 .OR. beta <0.0)THEN
        PRINT*, 'Error: A probability should be in the interval [0,1]'
        STOP 
      END IF

      DO q = 1, noise_data
        READ (10,*) channel_noise(q)
      
      ! Allocate arrays  
      ALLOCATE(weights(nodes, nodes), c1(nodes))
      ALLOCATE(S(nodes,Nstep), S1(nodes,Nstep), V(nodes,Nstep))
      ALLOCATE(m(nodes,Nstep), h(nodes,Nstep), n(nodes,Nstep), u(nodes,Nstep))
      ALLOCATE(W1(nodes,Nstep), dW1(nodes), dS1(nodes))
      ALLOCATE(SUM_ISI1(nodes), SUM_ISI_sq1(nodes)) 
      ALLOCATE(MEAN_ISI1(nodes), MEAN_ISI_sq1(nodes))
      ALLOCATE(sum_network_V(Nstep), sum_network2_V(Nstep), std_dev(Nstep))
  
      DO iseed = 1, realizations !Seed allocation for results reproducibility
          
        j = 1
        CALL RANDOM_SEED(SIZE = j)
        ALLOCATE(SEED(j))
        SEED = 12345 + (iseed-1)*1000
        CALL RANDOM_SEED(PUT = SEED)
        !PRINT*, "Current seed:", SEED(j)
        DEALLOCATE(SEED)
        !PRINT*, "Current seed:", SEED
      
    
      weights = 0  ! Initializing the weight matrix
      c1 = 0; S1 = 0.0d0 
      
      call cpu_time(tic)
      
      !Initial conditions
      DO k = 1, nodes
        CALL RANDOM_NUMBER(a11)
        CALL RANDOM_NUMBER(a22)
        CALL RANDOM_NUMBER(a33)
        CALL RANDOM_NUMBER(a44)
        V(k,1) = a11*(40 - 75.0d0); m(k,1) = a22*(1.0d0 - 0.0d0)
        h(k,1) = a33*(1.0d0 - 0.0d0); n(k,1) = a44*(1.0d0 - 0.0d0)
        u(k,1) = 0.01d0
      
        CALL RANDOM_NUMBER(a1)
        CALL RANDOM_NUMBER(a2)
        W1(k, 1) = mean_W1 + ( sd_W1 * (SQRT(-4.0d0*dt*DLOG(a1))*DCOS(2*pi*a2)) )
        !WRITE(14,*) W1(k,1)
      
        !Initializing the regular network topology
        DO j = k-neigh, k+neigh !Loop over the neighbours
          IF (k .NE. j) THEN ! Forbidding self-rewiring
            iv = j
            IF (j .LE. 0) iv = j + nodes ! Enforcing periodic boundary conditions
            IF (j .GT. nodes) iv = MOD(j,nodes) ! Enforcing periodic boundary conditions
      
            weights(k, iv) = 1; weights(iv, k) = 1 ! Making sure the weight matrix is diagonal
          END  IF
       END DO 
      
      END DO
      
      !Initializing the small-world topology
      DO k = 1, nodes
        DO j = k+1, k+neigh ! Sum over neighbours to the right only to avoid double-counting
      
          CALL RANDOM_NUMBER(rnd) ! A rewiring proposal is made for each connection of the node i
      
          IF (rnd .LT. beta) THEN !Rewiring is made with probability beta
            iv = j
            IF (j .GT. nodes) iv = MOD(j,nodes)
      
            weights(k,iv) = 0; weights(iv, k) = 0 ! Unwiring current connection 
            weights(k,k) = 2 ! Making sure diagonal elements are not chosen in the function PACK
            unwired_connections = PACK([(ind, ind=1,nodes)], weights(k,:) == 0 ) ! Select nodes available for rewiring
            
            CALL RANDOM_NUMBER(rnd)
      
            inode = INT(rnd*size(unwired_connections) + 1.0d0) ! Choose probalilistical the indice of the node that will be rewired
            iiv = unwired_connections(inode) ! Select the node that will be rewired
            weights(k, iiv) = 1; weights(iiv, k) = 1 ! Perform the rewiring
            weights(k,k) = 0 ! Set back the diagonal term to zero
      
          END IF
        END DO
      
      END DO
  
      ! Saving the initial graph
      !DO k = 1, nodes
      !  WRITE(14,*) weights(k, :) 
      !  print*, weights(k, :)  
      !END DO
      !CLOSE(14 )
      
      call cpu_time(tac)
      
      !PRINT*, "The total number of connection should be:", neigh*nodes
      !PRINT*, "The total number of connection is:", 0.5*sum(weights)
      !PRINT*

      
      t = 0.0d0
      
      call cpu_time(tic)
     
      sum_network_V = 0.0d0
      sum_network2_V = 0.0d0
    
      DO p = 1, Nstep - 1  !Runge Kutta loop starts
      
          DO k = 1, nodes
            alpha_m = 0.1d0*(V(k,p) + 40.0d0)/(1.0d0-EXP(-(V(k,p) + 40.d0)/10.0d0))
            alpha_h = 0.07d0*EXP(-(V(k,p) + 65.0d0)/20.0d0)
            alpha_n = 0.01d0*(V(k,p) + 55.0d0)/(1.0d0-EXP(-(V(k,p) + 55.0d0)/10.0d0))
            beta_m  = 4.0d0*EXP(-(V(k,p) + 65.0d0)/18.0d0)
            beta_h  = 1.0d0/(1.0d0 + EXP(-(V(k,p) + 35.0d0)/10.0d0))
            beta_n  = 0.125d0*EXP(-(V(k,p) + 65.0d0)/80.0d0)
            
            CALL RANDOM_NUMBER(a1)
            CALL RANDOM_NUMBER(a2)
            noise_V =  0.0d0 !SQRT(-4.0d0*sigma*dt*DLOG(a1))*DCOS(2*pi*a2)
  
            CALL RANDOM_NUMBER(a3)
            CALL RANDOM_NUMBER(a4)
            noise_m = SQRT(-4.0d0*dt*(alpha_m*beta_m)/((rho_Na*channel_noise(q))*(alpha_m + beta_m))*DLOG(a3))*DSIN(2*pi*a4)
  
            CALL RANDOM_NUMBER(a5)
            CALL RANDOM_NUMBER(a6)
            noise_h = SQRT(-4.0d0*dt*(alpha_h*beta_h)/((rho_Na*channel_noise(q))*(alpha_h + beta_h))*DLOG(a5))*DCOS(2*pi*a6)
  
            CALL RANDOM_NUMBER(a7)
            CALL RANDOM_NUMBER(a8)
            noise_n = SQRT(-4.0d0*dt*(alpha_n*beta_n)/((rho_K*channel_noise(q))*(alpha_n + beta_n))*DLOG(a7))*DSIN(2*pi*a8)
      
            ntk = 0.0d0
  
            ! Loop evaluating the synaptic current
            DO j = k+1, k+nodes-1 ! Loop over all the nodes to account for long-range connections
              iv = j
              IF (j .GT. nodes) iv = MOD(j,nodes) 
              
              IF(weights(k, iv) .NE. 0)THEN
      
                dS1(k) = S1(k, c1(k)) - S1(iv, c1(iv)) ! implement nearest spiking time of neighboring neuron iv
                !IF ( (p  .GT. 250000) .AND. MOD(p, 50000)==0) THEN
                  !print*, "t:", p, "i:", k, "i neigh:", iv
                  !print*, "i current # of spikes", c1(k), "i spiking time", S1(k, c1(k)), "time difference", dS1(k)
                  !print*, "i neigh current # of spikes", c1(iv), "i neigh spiking time", S1(iv, c1(iv)), "current time", p*dt 
                  !print*
                !end if
                
                !*******************************************************************************************************************
                ! Switch learning rules here
                !*******************************************************************************************************************
                ! !Hebbian learning rule
                 IF (dS1(k) .GE. 0.0d0) THEN   
                  dW1(k) = A_1 * EXP(-dS1(k)/tau_1)
                  !gstar = max_W1 
                ELSE
                  dW1(k) = - A_2 * EXP(dS1(k)/tau_2)
                  !gstar = min_W1
                END IF
                !W1(k,p+1) = W1(k,p) + ((gstar - W1(k,p)) * learning_rate  *  ABS(dW1(k))) ! multiplicative stdp
                 W1(k,p+1) = W1(k,p) + (learning_rate  *  dW1(k))    ! additive stdp
                 IF ( W1(k,p+1) .LT. min_W1 )THEN
                   W1(k,p+1) = min_W1
                 ELSE IF ( W1(k,p+1) .GT. max_W1 )THEN
                  W1(k,p+1) = max_W1
                 END IF
                 
                ! !Anti-Hebbian learning rule
                ! IF (dS1(k).GT.0.0d0) THEN   
                !  dW1(k) = - B_1 * EXP(-dS1(k)/tau_1)
                !   !gstar = min_W1 
                ! ELSE
                !   dW1(k) = - B_2 * dS1(k)/tau_2 * EXP(dS1(k)/tau_2)
                !   !gstar = max_W1
                ! END IF
                ! !W1(k,p+1) = W1(k,p) + ((gstar - W1(k,p)) * learning_rate  *  ABS(dW1(k))) ! multiplicative stdp
                !  W1(k,p+1) = W1(k,p) + (learning_rate  *  dW1(k))    ! additive stdp
                !      IF ( W1(k,p+1) .LT. min_W1  )THEN
                !        W1(k,p+1) = min_W1
                !      ELSE IF ( W1(k,p+1) .GT. max_W1  )THEN
                !        W1(k,p+1) = max_W1
                !      END IF
                !*******************************************************************************************************************
  
      
                !cnt = 0.0d0
                !DO c2 = 1, c1(iv)
                !  test_time = p*dt - S1(iv,c2) - tau_l  
                !  IF( test_time .GE. 0.0d0) THEN
                !    cnt = cnt + ( 1.0d0/(tau_d -tau_r) * ( EXP(- test_time/tau_d ) - EXP(- test_time/tau_r ) ) )
                !  ELSE 
                !    cnt = cnt + 0.0d0
                !  END IF
                !END DO 
                 
                ntk = ntk + weights(k, iv) * W1(k, p) * u(iv, p) * ( V(k,p) - V_r ) 
                !PRINT*, "weights(k,iv):", weights(k,iv) 
              END IF
            END DO
              
            p1=dt*f0(V(k,p), m(k,p), h(k,p), n(k,p), u(k,p),           t)
            v1=dt*f1(V(k,p), m(k,p), h(k,p), n(k,p), u(k,p), ntk,      t) + noise_V
            m1=dt*f2(V(k,p), m(k,p), h(k,p), n(k,p), u(k,p),           t) + noise_m
            h1=dt*f3(V(k,p), m(k,p), h(k,p), n(k,p), u(k,p),           t) + noise_h
            n1=dt*f4(V(k,p), m(k,p), h(k,p), n(k,p), u(k,p),           t) + noise_n
            
            p2=dt*f0(V(k,p) + v1/2.0d0, m(k,p) + m1/2.0d0, h(k,p) + h1/2.0d0, n(k,p) + n1/2.0d0, u(k,p) + p1/2.0d0,      &
            t + dt/2.0d0)
            v2=dt*f1(V(k,p) + v1/2.0d0, m(k,p) + m1/2.0d0, h(k,p) + h1/2.0d0, n(k,p) + n1/2.0d0, u(k,p) + p1/2.0d0, ntk, &
             t + dt/2.0d0) + noise_V
            m2=dt*f2(V(k,p) + v1/2.0d0, m(k,p) + m1/2.0d0, h(k,p) + h1/2.0d0, n(k,p) + n1/2.0d0, u(k,p) + p1/2.0d0,      & 
            t + dt/2.0d0) + noise_m
            h2=dt*f3(V(k,p) + v1/2.0d0, m(k,p) + m1/2.0d0, h(k,p) + h1/2.0d0, n(k,p) + n1/2.0d0, u(k,p) + p1/2.0d0,      & 
            t + dt/2.0d0) + noise_h
            n2=dt*f4(V(k,p) + v1/2.0d0, m(k,p) + m1/2.0d0, h(k,p) + h1/2.0d0, n(k,p) + n1/2.0d0, u(k,p) + p1/2.0d0,      &
            t + dt/2.0d0) + noise_n
                
            p3=dt*f0(V(k,p) + v2/2.0d0, m(k,p) + m2/2.0d0, h(k,p) + h2/2.0d0, n(k,p) + n2/2.0d0, u(k,p) + p2/2.0d0,      &
            t + dt/2.0d0)
            v3=dt*f1(V(k,p) + v2/2.0d0, m(k,p) + m2/2.0d0, h(k,p) + h2/2.0d0, n(k,p) + n2/2.0d0, u(k,p) + p2/2.0d0, ntk, &
            t + dt/2.0d0) + noise_V
            m3=dt*f2(V(k,p) + v2/2.0d0, m(k,p) + m2/2.0d0, h(k,p) + h2/2.0d0, n(k,p) + n2/2.0d0, u(k,p) + p2/2.0d0,      &
            t + dt/2.0d0) + noise_m
            h3=dt*f3(V(k,p) + v2/2.0d0, m(k,p) + m2/2.0d0, h(k,p) + h2/2.0d0, n(k,p) + n2/2.0d0, u(k,p) + p2/2.0d0,      &
            t + dt/2.0d0) + noise_h
            n3=dt*f4(V(k,p) + v2/2.0d0, m(k,p) + m2/2.0d0, h(k,p) + h2/2.0d0, n(k,p) + n2/2.0d0, u(k,p) + p2/2.0d0,      &
            t + dt/2.0d0) + noise_n
             
            p4=dt*f0(V(k,p) + v3, m(k,p) + m3, h(k,p) + h3, n(k,p) + n3, u(k,p) + p3,           t + dt)
            v4=dt*f1(V(k,p) + v3, m(k,p) + m3, h(k,p) + h3, n(k,p) + n3, u(k,p) + p3, ntk,      t + dt) + noise_V
            m4=dt*f2(V(k,p) + v3, m(k,p) + m3, h(k,p) + h3, n(k,p) + n3, u(k,p) + p3,           t + dt) + noise_m
            h4=dt*f3(V(k,p) + v3, m(k,p) + m3, h(k,p) + h3, n(k,p) + n3, u(k,p) + p3,           t + dt) + noise_h
            n4=dt*f4(V(k,p) + v3, m(k,p) + m3, h(k,p) + h3, n(k,p) + n3, u(k,p) + p3,           t + dt) + noise_n
      
            u(k,p + 1) = u(k,p) +  (p1 + 2.0d0*(p2 + p3) + p4)/6.0d0
            V(k,p + 1) = V(k,p) +  (v1 + 2.0d0*(v2 + v3) + v4)/6.0d0
            m(k,p + 1) = m(k,p) +  (m1 + 2.0d0*(m2 + m3) + m4)/6.0d0
            h(k,p + 1) = h(k,p) +  (h1 + 2.0d0*(h2 + h3) + h4)/6.0d0
            n(k,p + 1) = n(k,p) +  (n1 + 2.0d0*(n2 + n3) + n4)/6.0d0

          END DO  !Network loop K ends
    
        DO k = 1, nodes
          !If there is a spike, calculate the spiking time of this spike for each neuron (k) at the current time (p)
          IF (V(k,p) .GT. thres .AND. V(k,p + 1) .LT. thres) THEN
            c1(k) = c1(k) + 1
            S1(k, c1(k)) = p*dt
          END IF
          
          sum_network_V(p) = sum_network_V(p) + V(p, k)**2
          sum_network2_V(p) = sum_network2_V(p) + V(p, k)
         
        END DO
        
        mean_network_V = sum_network_V(p)/DBLE(nodes)
        mean_network2_V = ( sum_network2_V(p)/DBLE(nodes) )**2
        std_dev(p) = SQRT((mean_network_V -  mean_network2_V)/DBLE(nodes-1))

      t = t + dt
      END DO  !Runge kutta loop ends
      
      sum_std_dev = SUM(std_dev) 
      mean_std_dev(iseed) = sum_std_dev/DBLE(Nstep) 
      
        call cpu_time(tac)
           !From here, we start calulating the synaptic weigths W1(k,Nstep) and its standard deviation       
         Sum_1(iseed) = SUM(W1(:, Nstep))/DBLE(nodes) 
         s_deviation(iseed) = SQRT( (SUM( (W1(:, Nstep))**2 )/DBLE(nodes) - (Sum_1(iseed))**2)/dble(nodes-1))
  
        !From here we start calulating the coefficient of variation (CV) and syn_error
        S = 0.0d0
      DO k = 1, nodes
          c = 0; p = trans_time 
          DO WHILE(p .LE. Nstep)
            IF(V(k,p) .GT. thres .AND. V(k,p + 15) .LT. thres)THEN
              c = c + 1
              S(k,c) = p*dt
            END IF
            p = p + 15
          END DO
          !WRITE(*,*)'Node and number of spikes:', k, c
          
          IF (ALLOCATED(ISI1)) DEALLOCATE(ISI1)
          ALLOCATE(ISI1(k,c))
          ISI1 = 0.0D0; SUM_ISI1 = 0.0D0; SUM_ISI_sq1 = 0.0D0
          r = 2
          DO WHILE(r .LE. c)
            ISI1(k,r) = S(k,r) - S(k,r-1)
            SUM_ISI1(k) = SUM_ISI1(k) + ISI1(k,r)
            SUM_ISI_sq1(k) = SUM_ISI_sq1(k) + ISI1(k,r)**2
            r = r + 1
          END DO
          
          MEAN_ISI1(k) = SUM_ISI1(k)/dble(r-2)
          MEAN_ISI_sq1(k) = SUM_ISI_sq1(k)/dble(r-2)
      END DO ! k loop ends here, the mean over each neuron is calculated. Now we calculate the mean over all neuron
        
        SUM_MEAN_ISI1 = 0.0D0; SUM_MEAN_ISI_sq1 = 0.0D0
        DO k = 1, nodes
          SUM_MEAN_ISI1 = SUM_MEAN_ISI1 + MEAN_ISI1(k)
          SUM_MEAN_ISI_sq1 = SUM_MEAN_ISI_sq1 + MEAN_ISI_sq1(k)
        END DO
        
        MEAN1 =  SUM_MEAN_ISI1/dble(nodes)
        MEAN_sq1 = SUM_MEAN_ISI_sq1/dble(nodes)
        CV(iseed) = SQRT(MEAN_sq1 - MEAN1**2)/MEAN1
      END DO ! End loop on iseed realizations
      
   !Now we calculate (1) the average synaptic weight and its average std dev, (2) the Average CV over realization
        Sum_2 = 0.0d0;  s_deviation_2 = 0.0d0; CV_sum = 0.0d0; CV_sum2 = 0.0d0; syn_index_sum =0.0d0;syn_index_sum2 =0.0d0
        DO  iseed = 1, realizations 
          Sum_2 = Sum_2 + Sum_1(iseed) 
          s_deviation_2 = s_deviation_2 + s_deviation(iseed)
          CV_sum = CV_sum + CV(iseed)
          CV_sum2 = CV_sum2 + ( CV(iseed) )**2
          syn_index_sum = syn_index_sum + mean_std_dev(iseed)
          syn_index_sum2 = syn_index_sum2 + ( mean_std_dev(iseed) )**2
        END DO
        
        avg_w = Sum_2/DBLE(realizations)
        avg_sd = s_deviation_2/DBLE(realizations)
        CV_ave = CV_sum/DBLE(realizations)
        CV_sd = SQRT( (CV_sum2/dble(realizations) - CV_ave**2 )/dble(realizations - 1) )  
        synronization = syn_index_sum/DBLE(realizations)
        synronization_sd = SQRT ( (syn_index_sum2/dble(realizations) - synronization**2)/dble(realizations - 1) )
    
        WRITE(11,*)  channel_noise(q), avg_w, avg_sd
        WRITE(12,*)  channel_noise(q), CV_ave, CV_sd
        WRITE(13,*)  channel_noise(q), synronization, synronization_sd
        
      DEALLOCATE(weights, c1)
      DEALLOCATE(S, S1, V)
      DEALLOCATE(m, h, n ,u)
      DEALLOCATE(W1, dW1, dS1)
      DEALLOCATE(SUM_ISI1, SUM_ISI_sq1) 
      DEALLOCATE(MEAN_ISI1, MEAN_ISI_sq1)
      DEALLOCATE(sum_network_V, sum_network2_V, std_dev)

      END DO

      call cpu_time(finish)
      PRINT '("Total simulation time is: ",f25.5," seconds.")',finish-start
  
      
      CONTAINS
  
      DOUBLE PRECISION FUNCTION f0(x, y, w, z, u2, t)
          IMPLICIT NONE
          DOUBLE PRECISION, INTENT(in) :: x, y, w, z, u2, t
          f0=  5.0d0*(1.0d0 - u2)/( 1 + EXP( - (x + 3.0d0)/8.0d0  ) ) - u2
        END FUNCTION f0
      
        DOUBLE PRECISION FUNCTION f1(x, y, w, z, u2, syn, t)
          IMPLICIT NONE
          DOUBLE PRECISION, INTENT(in) :: x, y, w, z, u2, syn, t
          f1=(-g_Na*y**3*w*(x-V_Na)-g_K*z**4*(x-V_K)-g_L*(x-V_L) + I_0 - syn)/C_m    
        END FUNCTION f1
      
        DOUBLE PRECISION FUNCTION f2(x, y, w, z, u2, t)
          IMPLICIT NONE
          DOUBLE PRECISION, INTENT(in) :: x, y, w, z, u2, t
          f2=0.1d0*(x + 40.0d0)/(1.0d0-EXP(-(x + 40.d0)/10.0d0))*(1.0d0-y) &
          - 4.0d0*EXP(-(x + 65.0d0)/18.0d0)*y
        END FUNCTION f2
      
        DOUBLE PRECISION FUNCTION f3(x, y, w, z, u2, t)
          IMPLICIT NONE
          DOUBLE PRECISION, INTENT(in) :: x, y, w, z, u2, t
          f3=0.07d0*EXP(-(x + 65.0d0)/20.0d0)*(1.0d0-w) &
          - 1.0d0/(1.0d0 + EXP(-(x + 35.0d0)/10.0d0))*w
        END FUNCTION f3
      
        DOUBLE PRECISION FUNCTION f4(x, y, w, z, u2, t)
          IMPLICIT NONE
          DOUBLE PRECISION, INTENT(in) :: x, y, w, z, u2, t
          f4=0.01d0*(x + 55.0d0)/(1.0d0-EXP(-(x + 55.0d0)/10.0d0))*(1.0d0-z) &
          - 0.125d0*EXP(-(x + 65.0d0)/80.0d0)*z
        END FUNCTION f4
      
      END PROGRAM Network