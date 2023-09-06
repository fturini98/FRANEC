subroutine Cattura_DM()

    

    !PER VEDERE LE FUNZIONI DEI MODULI GUARDA IL FILE moduli.f90
  
    !Il modulo fisica serve a definire la matrice fisica G dove:
    !     G(1,k) e'      r/10^10          in cgs
    !     G(2,k) e' la   l/10^32          in cgs
    !     G(3,k) e' la  (P_totale)/10^17  in cgs
    !     G(4,k) e' la  T/10^6            in cgs     
    !     G(5,k) e' la  m/10^33           in cgs
    use fisica
  
    !Il modulo chimic sereve a definire la matrice della chimica XXX
    use chimic
  
    use chim
    !Il modulo chim serve per definire l' array XX, il passo temporale HT1 etc..
  
    implicit none
     
    integer :: k, counter=0
  
    !real,dimension(MELE,LIM) :: XXX
    
    write(*,*)"Il limite è:",LIM
    !write(*,*)"Piango"
    !write(*,*)"Il massimo del raggio è:",G(1,LIM)
    write(*,*)"Cambio 1" 
    do k = 1, LIM
        if (G(1,k)>0) then
            write(*,*)G(1,k)
            counter= counter +1
        endif
    end do
    write(*,*)"Speriamo Funzioni, Il counter è:",counter
      
  end subroutine Cattura_DM