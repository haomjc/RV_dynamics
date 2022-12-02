import numpy as np
from fh import fh
def rv_dyna(q_input,t,omega_hx, T_all):

    #均采用国际单位si
    
    P= 3000  #输入功率    w
    T_2 = 81.6   #输出扭矩： 应是N·m

    D_z = 400*10**(-3)
    Z_pie = 20
    
    Z = 20
    L = 20*10**(-3)
    
    b  = 5*10**(-3)          #渐开线齿宽
    L_pie = 1*10**(-3)
    n = 30     #rpm
    m = 1.5*10**(-3)
    z_1 = 15
    z_2  = 40
    z_g = 39
    z_b = z_g+1

    K_1 = 0.75
    #i = z_b*z_2/z_1+1
    i_16_5 = z_b*z_2/z_1+1
    i_15_6 = z_b*z_2/z_1
    i_65_1 = i_15_6/i_16_5    
    B = 10*10**(-3)
    d_z = 5*10**(-3)
    R_z = D_z/2
    G = 79380*10**6    #切变模量(pa)
    E = 212000*10**6  #弹性模量(Pa)
    #E = 2.07*10**7      #需要再查一下
    miu = 0.31 #泊松比
    
    r_1 = m*z_1/2  #齿轮
    r_2 = m*z_2/2
    
    e_c = K_1*R_z/z_b
    
    #===========输入轴的扭转刚度================
    
    L1 = 10+b/2+L_pie/2
    L2 = 3/2*B+L_pie/2
    L3 = 1/2*(B+L_pie)
    
    T_1= 0.78
    #T_1= 9.550*P/n   #单位是Nm???
    
    d1 = m*z_1
    L_inputShaft = L1+L2+L3     #10mm为输出端盘预留空间
    d_v = (L_inputShaft/(L1/d1**4+L2/(d1+5)**4+L3/(d1+10)**4))**(1/4)  #输入轴等效直径
    
    I_p = np.pi*d_v**4/32
    
    theta_si = T_1*(L_inputShaft)/(G*I_p)
    
    
    #k_01 = T_1/theta_si  #
    k_01 = 7.6769*10**4


    

    
    #=============渐开线啮合刚度================
    alpha = np.pi/9
    z_b = z_g+1
    h_a1 = 1*m
    h_a2 = 1*m
    d_a1 = m*z_1+2*h_a1
    d_a2 = m*z_2+2*h_a2
    alpha_a1 = np.arccos(m*z_1*np.cos(alpha)/d_a1)    
    alpha_a2 = np.arccos(m*z_2*np.cos(alpha)/d_a2) 
    epsilon_1 = z_1/(2*np.pi)*(np.tan(alpha_a1)-np.tan(alpha))
    epsilon_2 = z_2/(2*np.pi)*(np.tan(alpha_a2)-np.tan(alpha)) 
    
    n_p= 3  #行星轮数量
    
    F_t = 2*(9550*P/n)/(n_p*m*z_1)   #切向力 ？
    
    
    q = 0.04723+0.15551/z_1+0.25791/z_2
    c_th_pie = 1/q
    C_M = 0.8
    C_R = 1
    h_f  = 1.25*m
    C_B = (1+0.5*(1.2-h_f/m))
    epsilon_alpha = epsilon_1+epsilon_2
    c_pie = c_th_pie*C_M*C_R*C_B
    
    
    #c_r = (0.75*epsilon_alpha+0.25)*c_pie*1000  #总刚度,  N/m化为N/mm
    k_h= (0.75*epsilon_alpha+0.25)*c_pie
    #k_12  = c_r*b*10**6*10**3 #总刚度,  N/m化为N/mm? -> N/m
    k_hm = k_h*b*10**9

    #r_s = m*z_1/2
    #theta_sp = F_t/(c_r*r_s)

    #=================转臂轴承刚度======================   （存在过约束考虑并联）
    
    K_y = 2/np.pi*(1/K_1+(K_1**2-1)/(2*K_1**2)*np.log((1+K_1)/(1-K_1)))
    R_t = T_2*z_2*z_b/(6*i_16_5*z_1*K_1*z_b)
    R_x = T_2*z_2*z_b/(3*i_16_5*z_1*K_1*D_z)
    R_y = K_y*T_2*z_2*z_b/(3*i_16_5*z_1*K_1*D_z)
    R_abs = abs(R_t)+(R_x**2+R_y**2)**(1/2)
    
    #k_34 = 0.340*10**4*(R_abs**0.1*Z**0.9*L**0.8)*10**3   #
    k_34 = 3.7*10**8
    #delta_cb = R_abs/K_r
    a_sp = m*(z_1+z_2)/2
    #tau = z_g*delta_cb*i_65_1/(z_b*a_sp)
    
    #=================支承轴承刚度======================        
    
    R_pie = T_2/(2*n_p*(a_sp))   #受力存疑
    #k_54 = 0.340*10**4*(R_pie**0.1*Z_pie**0.9*L_pie**0.8)*10**3   #
    k_54 = 3.3*10**8
    #delta_sb = R_pie/K_r_pie
    #zeta = delta_sb/a_sp  

    #摆线轮刚度    ===================刚度过大====================
    
    phi_i = np.pi/2   #摆线轮当前转角
    
    _lambda=0.65
    #计算c_bz
    
    S = 1+K_1**2-2*K_1*np.cos(phi_i)
    T = K_1*(2+z_g)*np.cos(phi_i)-(1+(z_g+1)*K_1**2)
    #c_st = np.pi*B*E*D_z*S**(3/2)/(4*(1-miu**2)*(D_z*S**(3/2)+d_z*T+d_z*abs(T)))
    #c_st = np.pi*B*E*D_z*S**(3/2)/(4*(1-miu**2)*(D_z*S**(3/2)+2*d_z*abs(T)))  
    k_ni = np.pi*B*E*D_z*S**(3/2)/(4*(1-miu**2)*(D_z*S**(3/2)+2*d_z*abs(T))) 
    
    e = K_1*D_z/(2*z_b)
    L_i_pie = e*z_g*np.sin(phi_i)/(S**(1/2))
    #c_bz = _lambda*c_st*L_i_pie**2*(z_b/2)   
    #k_36 = _lambda*k_ni*L_i_pie**2*(z_b/2)   
    k_36 = 1.5*10**7
        #c_bz = c_bz+_lambda*c_st*(L_i_pie**2)*(z_b/2)
    #c_bz = c_bz/1000    #求平均扭转刚度
    #beta_H = T_g/c_bz
    #beta = i_65_1*z_b/z_g*beta_H    



    J_0 = 73.177*10**(-6)
    
    #ha=1;
    m1=np.pi*(m*z_1/2)**2*b*7800;
    m2=np.pi*(m*z_2/2)**2*b*7800;#单位kg    
    
    J_1=0.5*m1*(m*z_1/2)**2;  #单位kg*m**2
    J_2=0.5*m2*(m*z_2/2)**2;

    J_3 = 723.3358*10**(-6)   #摆线轮转动惯量
    J_5 = 896.7418*10**(-6)       #行星架转动惯量
    J_6 = 0.0200*10**(-6)    #滚针转动惯量   
    '''
    J_0 = 73.177*10**(-6)*10**(-6)
    J_1 = 2.7939*10**(-6)*10**(-6)    
    J_2 = 56.6449*10**(-6)*10**(-6)   #kg.m^2    
    J_3 = 723.3358*10**(-6)*10**(-6)    #摆线轮转动惯量
    J_5 = 896.7418*10**(-6)*10**(-6)        #行星架转动惯量
    J_6 = 0.0200*10**(-6)*10**(-6)         #滚针转动惯量

    J_0 = 73.177
    J_1 = 2.7939
    J_2 = 56.6449 #kg.m^2    

    J_3 = 723.3358   #摆线轮转动惯量
    J_5 = 896.7418       #行星架转动惯量
    J_6 = 0.0200    #滚针转动惯量       
    '''
    #===========输入轴的扭转阻尼================

    Xi = 0.005 #扭转阻尼比    
    c_01 = 2*Xi*np.sqrt(k_01/(1/J_0+1/J_1))
    #c_01 = 0.0045*10**(-3)
    #=============渐开线啮合阻尼================
    #几何计算
    r1=m*z_1/2;#分度圆半径
    r2=m*z_2/2;
    rb1=r1*np.cos(alpha);#基圆半径
    rb2=r2*np.cos(alpha);
    
    Xi_12 = 0.1  #渐开线阻尼比
    c_12 = 2*Xi_12*np.sqrt(k_hm*(r1**2*r2**2*J_1*J_2/(r1**2*J_1+r2**2*J_2)))
    #c_12 = 0.0214
    #=================摆线轮阻尼===============
    Xi_36 = 0.1    #齿轮材料啮合阻尼比
    #c_ni  = 2*Xi_36*np.sqrt(k_ni*())
    
    c_36= 4.8994*10**(-4)      #=========粗略值，后续需要改进


    c_34 = 0.001         #摆线轮 转臂轴承间阻尼
    c_54 = 0.001         #行星架与支撑轴承阻尼    

    #============时间标称尺度=================



    me=J_1*J_2/(J_1*rb2**2+J_2*rb1**2);#齿轮副等效质量
    #me=2/(1/J_1+1/J_2)
    
    #omega_d = np.sqrt(k_hm/me);      #齿轮固有频率，用于确定时间标称尺度
    #omega_d = 1
    #omega_n = omega_d
    omega_d = 4.7671*10**4

    #omega_h=omega_hx*omega_n;  #齿轮啮合频率
    #omega_h=omega_hx
    #omega_d = 1
    epsilon=0.1
    tau = t*omega_d
    Omega = omega_hx/omega_d
    #tau = t    
    #k_12 =  k_hm+epsilon*k_hm*np.cos(tau);   #时变刚度
    #k_12 =  k_hm+epsilon*k_hm*np.cos(z_1*Omega*tau);   #时变刚度,自己改的,啮合频率
    #k_12 =  k_hm
    #b_d = 10**(-6)*10**(4)
    b_d = 10**(-4)
    
    #=============无量纲化===============？？？需要改
    #m_0 = 0.54;m_1 = 0.0633;m_2 = 0.2786;m_3 = 0.5631;m_6 = 0.0045;m_5 = 4.9819    
    
    #========动力学方程表达式=================
    q =  q_input[:6]   #输入初始值
    dq =  q_input[6:]    #输入初始值
    print(q)
    ddq = np.zeros(6)
    
    x_01 = q[0]*b_d 
    x_12 = q[1]*b_d  
    x_12_pie = q[2]*b_d      
    x_23 = q[3]*b_d
    #x_23_pie = x[4]      
    #x_2pie_3 = x[5]
    x_2pie_3pie = q[4]*b_d 
    #x_36 = x[3]
    #x_3pie_6 = x[4]
    x_25 = q[5]*b_d
    #x_2pie_5pie = x[10]

    #============
    dx_01 = dq[0]*(b_d*omega_d) 
    dx_12 = dq[1]*(b_d*omega_d)   
    dx_12_pie = dq[2]*(b_d*omega_d)       
    dx_23 = dq[3]*(b_d*omega_d) 
    #dx_23_pie = dx[4]      
    #dx_2pie_3 = dx[5]
    dx_2pie_3pie = dq[4]*(b_d*omega_d)  
    #dx_36 = dx[3]
    #dx_3pie_6 = dx[4]
    dx_25 = dq[5]*(b_d*omega_d) 
    #dx_2pie_5pie = dx[10]
    #=============




    r_c_pie = R_z       #摆线轮节圆半径
    
    alpha_c = 9.5311*np.pi/180    #rad  摆线轮与针轮间啮合角
    
    
    zeta_12 = 10**(-6)   #综合误差
    zeta_12_pie = 10**(-6)
    
    #tau = t
    

    
    theta_0 = Omega*tau
    #theta_0 = omega_hx*t

    theta_c = theta_0/(z_2/z_1)
    #theta_c = theta_0/2.6  #曲柄轴自转角度    
    #theta_g = theta_0/105   #摆线轮自转角度
    theta_g = theta_0/i_16_5
    theta_1 = theta_0 - x_01
    
    theta_2 = (r_1*theta_1-x_12)/r_2
    theta_2_pie = (r_1*theta_1-x_12_pie)/r_2
    
    theta_3 = (theta_2 - x_23)/(i_16_5/(z_2/z_1))
    theta_3_pie = (theta_2_pie - x_2pie_3pie)/(i_16_5/(z_2/z_1))
    #theta_3 = x_36
    #theta_3_pie = x_3pie_6
    
    theta_5 = (theta_2 - x_25)/(i_16_5/(z_2/z_1))
    
    theta_3c = r_2/r_1*theta_2/(z_2/z_1) #自己推导的
    #theta_3c = theta_1/2.6
    #theta_3c = r_2/r_1*theta_2/2.6    
    
    theta_3c_pie = r_2/r_1*theta_2_pie/(z_2/z_1)
    #theta_3c_pie = r_2/r_1*theta_2_pie/2.6    
    #theta_3c_pie = theta_1/2.6
    #==========================
    d_theta_0 = omega_hx    #与theta_0要一致
    #d_theta_0 = Omega
    d_theta_c = d_theta_0/(z_2/z_1)  #曲柄轴自转角度
    #d_theta_c = d_theta_0/2.6  #曲柄轴自转角度    
    d_theta_g = d_theta_0/i_16_5   #摆线轮自转角度

    d_theta_1 = d_theta_0 - dx_01
    
    d_theta_2 = (r_1*d_theta_1-dx_12)/r_2
    d_theta_2_pie = (r_1*d_theta_1-dx_12_pie)/r_2
    
    d_theta_3 = (d_theta_2 - dx_23)/(i_16_5/(z_2/z_1))
    d_theta_3_pie = (d_theta_2_pie - dx_2pie_3pie)/(i_16_5/(z_2/z_1))

    #d_theta_3 = dx_36
    #d_theta_3_pie = dx_3pie_6
    d_theta_5 = (d_theta_2 - dx_25)/(i_16_5/(z_2/z_1))
    
    d_theta_3c = r_2/r_1*d_theta_2/(z_2/z_1)  #自己推导的
    #d_theta_3c = r_2/r_1*d_theta_2/2.6  #自己推导的    
    d_theta_3c_pie = r_2/r_1*d_theta_2_pie/(z_2/z_1)  
    #d_theta_3c_pie = r_2/r_1*d_theta_2_pie/2.6  
    
    k_12 =  k_hm+epsilon*k_hm*np.cos(z_1*theta_1);   #时变刚度,自己改的,啮合频率    
    
    
    f_12 = (k_12*(fh(x_12)-zeta_12)+c_12*dx_12)
    
    f_12_pie = (k_12*(fh(x_12_pie)-zeta_12_pie)+c_12*dx_12_pie)
    #======================
    i = np.array([0, 0, 1, 1]);j = np.array([0, 1, 0, 1]);theta_3cj = np.array([theta_3c,theta_3c_pie, theta_3c,theta_3c_pie ]);
    theta_3j =np.array([ theta_3, theta_3_pie, theta_3, theta_3_pie]);theta_2i = np.array([theta_2, theta_2, theta_2_pie, theta_2_pie])
    d_theta_3cj = np.array([d_theta_3c, d_theta_3c_pie, d_theta_3c, d_theta_3c_pie]);d_theta_3j =np.array([d_theta_3 , d_theta_3_pie, d_theta_3 , d_theta_3_pie]);
    d_theta_2i =np.array([d_theta_2, d_theta_2, d_theta_2_pie, d_theta_2_pie])
    
    tau_3j_4iY = e_c*(theta_3cj-theta_c)*np.sin(theta_c+np.pi*j)+a_sp*(theta_3j-theta_g)*np.sin(theta_g+np.pi*i)
    zeta_3j_4iY = 10**(-6)   #综合误差  
    sigma_2i_3jY = e_c*(theta_2i-theta_c)*np.sin(theta_c+np.pi*j)
    delta_3j_4iY = sigma_2i_3jY-tau_3j_4iY-zeta_3j_4iY
    
    d_sigma_2i_3jY = e_c*(d_theta_2i-d_theta_c)*np.sin(theta_c+np.pi*j)+e_c*(theta_2i-theta_c)*np.cos(theta_c+np.pi*j)*d_theta_c
    d_tau_3j_4i_Y = e_c*(d_theta_3cj-d_theta_3c)*np.sin(theta_c+np.pi*j)+e_c*(theta_3cj-theta_c)*np.cos(theta_c+np.pi*j)*d_theta_c+ \
    a_sp*(d_theta_3j-d_theta_g)*np.sin(theta_g+np.pi*i)+a_sp*(theta_3j-theta_g)*np.cos(theta_g+np.pi*i)*d_theta_g
    d_delta_3j_4iY = d_sigma_2i_3jY - d_tau_3j_4i_Y
    #d_delta_3j_4iY = e_c*(theta_2i-theta_c)*np.cos(theta_c+np.pi*j)+  e_c*(theta_3cj-theta_c)*np.cos(theta_c+np.pi*j)+a_sp*(theta_3j-theta_g)*np.cos(theta_g+np.pi*i)
    #F_34Y, F_3pie_4Y, F_3_4pieY, F_3pie_4pieY = k_34*delta_3j_4iY+c_34*d_delta_3j_4iY
    F_34Y, F_3pie_4Y, F_3_4pieY, F_3pie_4pieY = (k_34*delta_3j_4iY+c_34*d_delta_3j_4iY)
    #print(F_34Y, F_3pie_4Y, F_3_4pieY, F_3pie_4pieY)
    #====================X方向===================
    tau_3j_4iX = e_c*(theta_3cj-theta_c)*np.cos(theta_c+np.pi*j)-a_sp*(theta_3j-theta_g)*np.cos(theta_g+np.pi*i)
    zeta_3j_4iX = 10**(-6)   #综合误差 
    sigma_2i_3jX = e_c*(theta_2i-theta_c)*np.cos(theta_c+np.pi*j)
    delta_3j_4iX = sigma_2i_3jX-tau_3j_4iX-zeta_3j_4iX
    
    
    d_tau_3j_4i_X = e_c*(d_theta_3cj-d_theta_3c)*np.cos(theta_c+np.pi*j)-e_c*(theta_3cj-theta_c)*np.sin(theta_c+np.pi*j)*d_theta_c- \
    a_sp*(d_theta_3j-d_theta_g)*np.cos(theta_g+np.pi*i)+a_sp*(theta_3j-theta_g)*np.sin(theta_g+np.pi*i)*d_theta_g
    d_sigma_2i_3jX = e_c*(d_theta_2i-d_theta_c)*np.cos(theta_c+np.pi*j)-e_c*(theta_2i-theta_c)*np.sin(theta_c+np.pi*j)*d_theta_c    
    d_delta_3j_4iX = d_sigma_2i_3jX - d_tau_3j_4i_X
    #d_delta_3j_4iX = -e_c*(theta_2i-theta_c)*np.sin(theta_c+np.pi*j)-e_c*(theta_3cj-theta_c)*np.sin(theta_c+np.pi*j)+a_sp*(theta_3j-theta_g)*np.sin(theta_g+np.pi*i)
    F_34X, F_3pie_4X, F_3_4pieX,  F_3pie_4pieX= (k_34*delta_3j_4iX+c_34*d_delta_3j_4iX)
    #F_34X, F_3pie_4X, F_3_4pieX,  F_3pie_4pieX= k_34*delta_3j_4iX+c_34*d_delta_3j_4iX
    #==========================F_63==========================
    
    j =np.array([0, 1]);theta_3j =np.array([theta_3, theta_3_pie]) ;theta_3cj =np.array([theta_3c, theta_3c_pie]) ;d_theta_3j =np.array([d_theta_3 , d_theta_3_pie]);d_theta_3cj = np.array([d_theta_3c, d_theta_3c_pie])
    alpha_k = np.pi/2-alpha_c+np.pi*j        #？？？ 摆线轮与针轮间法向啮合角
    tau_3j_6 = -r_c_pie*(theta_3j-theta_g)*np.sin(alpha_k)+e_c*(theta_3cj-theta_c)*np.sin(alpha_k)
    zeta_3j_6 = 10**(-6)
    delta_3j_6 = tau_3j_6 - zeta_3j_6
    d_delta_3j_6 = -r_c_pie*(d_theta_3j-d_theta_g)*np.sin(alpha_k)+e_c*(d_theta_3cj-d_theta_c)*np.sin(alpha_k)
    F_63, F_63_pie = (k_36*delta_3j_6 + c_36*d_delta_3j_6)
    #F_63, F_63_pie = k_36*delta_3j_6 + c_36*d_delta_3j_6    
    #======================F_54====================
    zeta_4i5_Y = 10**(-6)
    zeta_4i5_X = 10**(-6)    
    i = np.array([0, 1]);
    tau_54i_Y = a_sp*(theta_5-theta_g)*np.sin(theta_g+np.pi*i)
    delta_54i_Y = -tau_54i_Y-zeta_4i5_Y
    d_delta_54i_Y = -a_sp*(d_theta_5-d_theta_g)*np.sin(theta_g+np.pi*i)-a_sp*(theta_5-theta_g)*np.cos(theta_g+np.pi*i)*d_theta_g
    F_54_Y , F_54pie_Y= (k_54*delta_54i_Y+c_54*d_delta_54i_Y)
    #F_54_Y , F_54pie_Y= k_54*delta_54i_Y+c_54*d_delta_54i_Y    
    #========X方向=================   
    tau_54i_X = -a_sp*(theta_5-theta_g)*np.cos(theta_g+np.pi*i)
    delta_54i_X = -tau_54i_X-zeta_4i5_X
    d_delta_54i_X = a_sp*(d_theta_5-d_theta_g)*np.cos(theta_g+np.pi*i)-a_sp*(theta_5-theta_g)*np.sin(theta_g+np.pi*i)*d_theta_g
    F_54_X, F_54pie_X = (k_54*delta_54i_X+c_54*d_delta_54i_X)
    #print(F_54_X,F_54pie_X )

    #F_54_X, F_54pie_X = k_54*delta_54i_X+c_54*d_delta_54i_X

    #=================
                
    ddx_01 = (T_1 - k_01*x_01-c_01*dx_01-(J_0/J_1)*(k_01*x_01+c_01*dx_01-r_1*f_12-r_1*f_12_pie))/J_0
    
    ddx_12 = (k_01*x_01+c_01*dx_01-r_1*f_12-r_1*f_12_pie-(J_1*r_2/(J_2*r_1))*(r_2*f_12-(F_34Y-F_3pie_4Y)*e_c*np.sin(theta_c)-(F_34X-F_3pie_4X)*e_c*np.cos(theta_c)))/(J_1/r_1)
    
    ddx_12_pie = (k_01*x_01+c_01*dx_01-r_1*f_12-r_1*f_12_pie-(J_1*r_2/(J_2*r_1))*(r_2*f_12_pie-(F_3_4pieY-F_3pie_4pieY)*e_c*np.sin(theta_c)-(F_3_4pieX-F_3pie_4pieX)*e_c*np.cos(theta_c)))/(J_1/r_1)
    
    ddx_23 = (r_2*f_12-(F_34Y-F_3pie_4Y)*e_c*np.sin(theta_c)-(F_34X-F_3pie_4X)*e_c*np.cos(theta_c)-i_16_5/(z_2/z_1)*(J_2/J_3)*(F_63*r_c_pie*np.cos(alpha_c)-F_34X*a_sp*np.cos(theta_g)+F_3_4pieX*a_sp*np.cos(theta_g)+F_34Y*a_sp*np.sin(theta_g)-F_3_4pieY*a_sp*np.sin(theta_g)))/J_2
    
    #ddx_23_pie = (r_2*f_12-(F_34Y-F_3pie_4Y)*e_c*np.sin(theta_c)-(F_34X-F_3pie_4X)*e_c*np.cos(theta_c)-(J_2/J_3)*(F_63_pie*r_c_pie*np.cos(alpha_c)-F_3pie_4X*a_sp*np.cos(theta_g)+F_3pie_4pieX*a_sp*np.cos(theta_g+F_3pie_4Y*a_sp*np.sin(theta_g)-F_3pie_4pieY*a_sp*np.sin(theta_g))))/J_2    
    
    #ddx_2pie_3 = (r_2*f_12_pie-(F_3_4pieY-F_3pie_4pieY)*e_c*np.sin(theta_c)-(F_3_4pieX-F_3pie_4pieX)*e_c*np.cos(theta_c)-(J_2/J_3)*(F_63*r_c_pie*np.cos(alpha_c)-F_34X*a_sp*np.cos(theta_g)+F_3_4pieX*a_sp*np.cos(theta_g+F_34Y*a_sp*np.sin(theta_g)-F_3_4pieY*a_sp*np.sin(theta_g))))/J_2     
    
    ddx_2pie_3pie = (r_2*f_12_pie-(F_3_4pieY-F_3pie_4pieY)*e_c*np.sin(theta_c)-(F_3_4pieX-F_3pie_4pieX)*e_c*np.cos(theta_c)-i_16_5/(z_2/z_1)*(J_2/J_3)*(F_63_pie*r_c_pie*np.cos(alpha_c)-F_3pie_4X*a_sp*np.cos(theta_g)+F_3pie_4pieX*a_sp*np.cos(theta_g)+F_3pie_4Y*a_sp*np.sin(theta_g)-F_3pie_4pieY*a_sp*np.sin(theta_g)))/J_2     
    
    #ddx_36 = (F_63*r_c_pie*np.cos(alpha_c)-F_34X*a_sp*np.cos(theta_g)+F_3_4pieX*a_sp*np.cos(theta_g)+F_34Y*a_sp*np.sin(theta_g)-F_3_4pieY*a_sp*np.sin(theta_g))/J_3
    
    #ddx_3pie_6 = (F_63_pie*r_c_pie*np.cos(alpha_c)-F_3pie_4X*a_sp*np.cos(theta_g)+F_3pie_4pieX*a_sp*np.cos(theta_g)+F_3pie_4Y*a_sp*np.sin(theta_g)-F_3pie_4pieY*a_sp*np.sin(theta_g))/J_3
    
    ddx_25 = (-(F_54pie_X-F_54_X)*a_sp*np.cos(theta_g)-(F_54_Y-F_54pie_Y)*a_sp*np.sin(theta_g)+T_2+i_16_5/(z_2/z_1)*(J_5/J_2)*(r_2*f_12-(F_34Y-F_3pie_4Y)*e_c*np.sin(theta_c)-(F_34X-F_3pie_4X)*e_c*np.cos(theta_c)))/J_5
    
    #ddx_2pie_5 = (-(F_54pie_X-F_54_X)*a_sp*np.cos(theta_g)-(F_54_Y-F_54pie_Y)*a_sp*np.sin(theta_g)+T_2+(J_5/J_2)*(r_2*f_12_pie-(F_3_4pieY-F_3pie_4pieY)*e_c*np.sin(theta_c)-(F_3_4pieX-F_3pie_4pieX)*e_c*np.cos(theta_c)))/J_5

    #===========test==========

    ddq[0]  =   ddx_01/(b_d*omega_d**2) 
    ddq[1]  =   ddx_12/(b_d*omega_d**2)  
    ddq[2]  =     ddx_12_pie/(b_d*omega_d**2) 
    ddq[3]  =   ddx_23/(b_d*omega_d**2) 
    #ddx[4]   =       ddx_23_pie 
    #ddx[5] =   ddx_2pie_3 
    ddq[4]   =  ddx_2pie_3pie/(b_d*omega_d**2)  
    #ddx[3]  =    ddx_36 
    #ddx[4]    =  ddx_3pie_6 
    ddq[5]  =   ddx_25/(b_d*omega_d**2) 
    #ddx[10]  =    ddx_2pie_5 
    print(str(t/T_all))

    return np.hstack([dq, ddq])
