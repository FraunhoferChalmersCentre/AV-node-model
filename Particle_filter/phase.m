function PHI=phase(G)

% Calculates the phase

%   L. Ljung 10-2-86
%   Copyright (c) 1986-98 by The MathWorks, Inc.
%   $Revision: 2.3 $  $Date: 1997/12/02 03:41:32 $

PHI=atan2(imag(G),real(G));
N=length(PHI);
DF=PHI(1:N-1)-PHI(2:N);
I=find(abs(DF)>3.5);
for i=I
    if i~=0
        PHI=PHI+2*pi*sign(DF(i))*[zeros(1,i) ones(1,N-i)];
    end

end

end