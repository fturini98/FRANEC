subroutine RESNUC(JFF)                                             
  use interfaccia
  use second
  use sistema
  use numer

  implicit none

  integer :: jff

  integer,parameter :: IPRQUI = 0
  real :: aa, ab, dd, ff, h
  real :: v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11, v12
  

  AA = C(1,1)/C(1,2)                                                  
  AB = C(1,4)/C(1,2)                                                  
  DD = C(2,4)/C(2,3)  
  FF = ((C(3,4)-AB*C(3,2))-DD*C(3,3))/(C(3,1)-AA*C(3,2))              
  H = (C(4,4)-AB*C(4,2)-DD*C(4,3))-FF*(C(4,1)-AA*C(4,2))  
  V1 = ((ALF1(4)-AB*ALF1(2)-DD*ALF1(3))-FF*(ALF1(1)-AA*ALF1(2)))/H    
  V2 = ((BET1(4)-AB*BET1(2)-DD*BET1(3))-FF*(BET1(1)-AA*BET1(2)))/H    
  V3 = ((GAM1(4)-AB*GAM1(2)-DD*GAM1(3))-FF*(GAM1(1)-AA*GAM1(2)))/H    
  H = C(3,1)-AA*C(3,2)                                                
  V4 = ((ALF1(1)-AA*ALF1(2))-V1*(C(4,1)-AA*C(4,2)))/H                 
  V5 = ((BET1(1)-AA*BET1(2))-V2*(C(4,1)-AA*C(4,2)))/H                 
  V6 = ((GAM1(1)-AA*GAM1(2))-V3*(C(4,1)-AA*C(4,2)))/H                 
  V7 = (ALF1(3)-C(3,3)*V4-C(4,3)*V1)/C(2,3)                           
  V8 = (BET1(3)-C(3,3)*V5-C(4,3)*V2)/C(2,3)                           
  V9 = (GAM1(3)-C(3,3)*V6-C(4,3)*V3)/C(2,3)                           
  V10 = (ALF1(2)-C(3,2)*V4-C(4,2)*V1)/C(1,2)                          
  V11  =(BET1(2)-C(3,2)*V5-C(4,2)*V2)/C(1,2)                          
  V12 = (GAM1(2)-C(3,2)*V6-C(4,2)*V3)/C(1,2)                          
  H = V3*V5-V6*V2                                                     
  ALF(3,JFF) = -(V3*V4-V6*V1)/H                                       
  BET(3,JFF) = V3/H                                                   
  GAM(3,JFF) = -V6/H                                                  
  ALF(4,JFF) = -(V1+V2*ALF(3,JFF))/V3                                 
  BET(4,JFF) = -V2*BET(3,JFF)/V3                                      
  GAM(4,JFF) = (1.-V2*GAM(3,JFF))/V3                                  
  ALF(1,JFF) = V10+V11*ALF(3,JFF)+V12*ALF(4,JFF)                      
  BET(1,JFF) = V11*BET(3,JFF)+V12*BET(4,JFF)                          
  GAM(1,JFF) = V11*GAM(3,JFF)+V12*GAM(4,JFF)                          
  ALF(2,JFF) = V7+V8*ALF(3,JFF)+V9*ALF(4,JFF)                         
  BET(2,JFF) = V8*BET(3,JFF)+V9*BET(4,JFF)                            
  GAM(2,JFF) = V8*GAM(3,JFF)+V9*GAM(4,JFF)                            

  if(IPRALL == 0 .and. IPRQUI == 0) return       
  
  write(2,111) ALF(1,JFF),BET(1,JFF),GAM(1,JFF),ALF(2,JFF),BET(2,JFF), &
       GAM(2,JFF),ALF(3,JFF),BET(3,JFF),GAM(3,JFF),ALF(4,JFF),BET(4,JFF), &
       GAM(4,JFF)                                                       
111 format(1X,1P,12E10.3)                                             
  return                                                            
end subroutine RESNUC
