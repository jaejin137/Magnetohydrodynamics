c     function to obtain system time
c
      function systime()

      double precision:: systime
      integer:: timearr(8)

      call date_and_time(values = timearr)
      systime = timearr(5)*3600+timearr(6)*60+timearr(7)+.001*timearr(8)
      
      end
