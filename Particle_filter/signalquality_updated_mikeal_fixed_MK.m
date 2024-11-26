function [Estimated_signal,Estimated_First_harmonic,Estimated_Second_harmonic,Estimated_Amp1,Estimated_Amp2,Estimated_freq,SQI,global_w_T]=...
    signalquality_updated_mikeal_fixed_MK(s3, lc1, lg1, delta_w, Fs)


% Mostafa Abdollahpur's code for paper DOI 10.3389/fphys.2022.976925

% Note, this code was written by Mostafa Abdollahpur and slightly changed
% by Mattias Karlsson.

% Lt= length of whole signal .
% lg1=(s) the window length for finding global frequency.
% lg2= (s)the overlap for window global frequency.
% lc1= (s)the window length for finding local frequency.
% lc2= (s)the overlap for window local frequency.

signal=hilbert(s3);
Lt=length(signal);


p1 = fix((length(signal))/(lg1*Fs));
base1=cell(1,p1);
base2=cell(1,p1);
ColumnMean_a_teta_a1=cell(1,p1);
ColumnMean_a_teta_a2=cell(1,p1);
base_freq=cell(1,p1);

ss4=cell(1,p1);

N=lg1*Fs-lc1*Fs+1; 

B=zeros(N,lg1*Fs);
B1=zeros(N,lg1*Fs);

Amplitude1=zeros(N,lg1*Fs);
Amplitude2=zeros(N,lg1*Fs);
Phase1=zeros(N,lg1*Fs);
Phase2=zeros(N,lg1*Fs);
pp1=linspace(0,lc1*Fs-1,lc1*Fs);

disp('Running Abdollahpur model for estimating frequency component in f-waves')
for k=1:p1
 
    disp([num2str(k), ' of ', num2str(p1)]) 

    signal1 = signal(1+(k-1)*(lg1*Fs):k*lg1*Fs,1);
     
    w0=@(w)norm((signal1-(Z_matrix(w,signal1)*inv(ctranspose(Z_matrix(w,signal1))*Z_matrix(w,signal1)))*(ctranspose(Z_matrix(w,signal1))*signal1)));
        
    w_g=[2*pi*4/Fs,2*pi*12/Fs];
    
    ifq=mean(instfreq(signal1,Fs));
    global_w=ifq;
    
    global_w_T(k)=global_w;
    i2=global_w-delta_w;
    i3=global_w+delta_w;

    for j=1:N
    
        signal2=signal1(j:(j-1)+lc1*Fs,1); 
        w0=@(w)norm((signal2-(Z_matrix(w,signal2)*inv(ctranspose(Z_matrix(w,signal2))*Z_matrix(w,signal2)))*(ctranspose(Z_matrix(w,signal2))*signal2)));
            
         w_l=[2*pi*i2/Fs,2*pi*i3/Fs];
            
        [asa_l,fval2] = fminbnd(w0,w_l(1),w_l(2));
        
        local_w=asa_l*Fs/(2*pi);
         
        B(j,j:j+lc1*Fs-1)=local_w;
        B1(j,j:j+lc1*Fs-1)=local_w*2*pi/Fs;
        pp2(j,j:j+lc1*Fs-1)=pp1;
        
        local_w1=2*pi*local_w/Fs;
        z1=Z_matrix(local_w1,signal2);
        z2=ctranspose(z1);
        a_teta=(inv(z2*z1))*z2*signal2;
        Amplitude1(j,j:j+lc1*Fs-1)=abs(a_teta(1,1));
        Amplitude2(j,j:j+lc1*Fs-1)=abs(a_teta(2,1));
        
        Phase1(j,j:j+lc1*Fs-1)=phase(a_teta(1,1)); 
        Phase2(j,j:j+lc1*Fs-1)=phase(a_teta(2,1));
        
        signal2=0;
end

    ColumnMean_Freq = sum(B,1) ./ sum(B~=0,1);
    ColumnMean_a_teta_a1{1,k} = sum(Amplitude1,1) ./ sum(Amplitude1~=0,1);
    ColumnMean_a_teta_a2{1,k}= sum(Amplitude2,1) ./ sum(Amplitude2~=0,1);
    
    cc1=B1.*pp2+Phase1;
    cc2=(2*B1).*pp2+Phase2;
    d1=find(cc1==0);
    d2=find(cc2==0);
    cc1(d1)=pi/2;
    cc2(d2)=pi/2;
    q1=cos(cc1);
    q2=cos(cc2);
    
    zz1=zeros(N,length(q1));
    zz2=zeros(N,length(q2));
    zz1(find(abs(q1)>=0.000001))=q1(find(abs(q1)>=0.000001));
    zz2(find(abs(q2)>=0.000001))=q2(find(abs(q2)>=0.000001));
    ColumnMean_cc1 = sum(zz1,1) ./ sum(zz1~=0,1);
    ColumnMean_cc2 = sum(zz2,1) ./ sum(zz2~=0,1); 
    
    
    base1{1,k}=ColumnMean_cc1;
    base2{1,k}=ColumnMean_cc2;
    base_freq{1,k}=ColumnMean_Freq;

end


%% CVX seperated

for z = 1:p1

    signal3 = s3(1+(z-1)*(lg1*Fs):z*lg1*Fs,1);  
    base_a1=cell2mat(base1(1,z));
    base_a2=cell2mat(base2(1,z));
    
    n=length(signal3);
    B=zeros(n,n-1);
    for i=2:n
        B(i-1,i-1)=1;
        B(i,i-1)=-1;
    end
    delta_z=zeros(1,length(signal3)-1);
    delta_z(1,:)=0.01*std(signal3);
    cvx_begin
        variable z1(n);
        minimize (norm(base_a1'.*z1-signal3));
        subject to
        B'*z1 <= delta_z';
    cvx_end
    
    cvx_begin
        variable z2(n);
        minimize (norm(base_a2'.*z2-signal3));
        subject to
        B'*z2 == 0;
    cvx_end
    
    ss1=base_a1'.*z1;
    ss2=base_a2'.*z2;
    
    ss3=ss1+ss2;
        
    e=signal3-ss3;
    
    sqi(1+(z-1)*(lg1*Fs):z*lg1*Fs,1)=1-(std(e)/std(signal3));
      
    ss4{1,z}= ss3;
    ss10{1,z}=ss1;
    ss20{1,z}=ss2;
    
    Estimated_Amp1{1,z}=z1;
    Estimated_Amp2{1,z}=z2;
end

Estimated_signal=ss4;
Estimated_First_harmonic=ss10;
Estimated_Second_harmonic=ss20;
Estimated_Amp1=Estimated_Amp1;
Estimated_Amp2=Estimated_Amp2;
Estimated_freq=base_freq;
SQI=sqi;

end
