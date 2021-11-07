% 读取格式为Format212的心电数据
% 标注也显示在图中，ANNOT中为标注，相应的时间在ATRTIME中。
% 标注存储为数字，其含义可在www.physionet.org的 "ecgcodes.h"中。
% -------------------------------------------------------------------------
clc; clear all;

%------ SPECIFY DATA ------------------------------------------------------
PATH= 'C:\Users\wang\Desktop\Pynq_temp\MIT-BIH Arrhythmia Database'; % path, where data are saved
HEADERFILE= '201.hea';      % header-file in text format
ATRFILE= '201.atr';         % attributes-file in binary format
DATAFILE='201.dat';         % data-file
SAMPLES2READ=30000;         % number of samples to be read
                            % in case of more than one signal:
                            % 2*SAMPLES2READ samples are read

%------ LOAD HEADER DATA --------------------------------------------------
fprintf(1,'\\n$> WORKING ON %s ...\n', HEADERFILE);
signalh= fullfile(PATH, HEADERFILE);
fid1=fopen(signalh,'r');
z= fgetl(fid1);
A= sscanf(z, '%*s %d %d %d',[1,3]);
nosig= A(1);  % number of signals
sfreq=A(2);   % sample rate of data
clear A;
for k=1:nosig
    z= fgetl(fid1);
    A= sscanf(z, '%*s %d %d %d %d %d',[1,5]);
    dformat(k)= A(1);           % format; here only 212 is allowed
    gain(k)= A(2);              % number of integers per mV
    bitres(k)= A(3);            % bitresolution
    zerovalue(k)= A(4);         % integer value of ECG zero point
    firstvalue(k)= A(5);        % first integer value of signal (to test for errors)
end;
fclose(fid1);
clear A;

%------ LOAD BINARY DATA --------------------------------------------------
if dformat~= [212,212], error('this script does not apply binary formats different to 212.'); end;
signald= fullfile(PATH, DATAFILE);            % data in format 212
fid2=fopen(signald,'r');
A= fread(fid2, [3, SAMPLES2READ], 'uint8')';  % matrix with 3 rows, each 8 bits long, = 2*12bit
fclose(fid2);
M2H= bitshift(A(:,2), -4);
M1H= bitand(A(:,2), 15);
PRL=bitshift(bitand(A(:,2),8),9);     % sign-bit
PRR=bitshift(bitand(A(:,2),128),5);   % sign-bit
M( : , 1)= bitshift(M1H,8)+ A(:,1)-PRL;
M( : , 2)= bitshift(M2H,8)+ A(:,3)-PRR;
if M(1,:) ~= firstvalue, error('inconsistency in the first bit values'); end;
switch nosig
case 2
    M( : , 1)= (M( : , 1)- zerovalue(1))/gain(1);
    M( : , 2)= (M( : , 2)- zerovalue(2))/gain(2);
    TIME=(0:(SAMPLES2READ-1))/sfreq;
case 1
    M( : , 1)= (M( : , 1)- zerovalue(1));
    M( : , 2)= (M( : , 2)- zerovalue(1));
    M=M';
    M(1)=[];
    sM=size(M);
    sM=sM(2)+1;
    M(sM)=0;
    M=M';
    M=M/gain(1);
    TIME=(0:2*(SAMPLES2READ)-1)/sfreq;
otherwise  % this case did not appear up to now!
    % here M has to be sorted!!!
    disp('Sorting algorithm for more than 2 signals not programmed yet!');
end;
clear A M1H M2H PRR PRL;
fprintf(1,'\\n$> LOADING DATA FINISHED \n');

%------ LOAD ATTRIBUTES DATA ----------------------------------------------
atrd= fullfile(PATH, ATRFILE);      % attribute file with annotation data
fid3=fopen(atrd,'r');
A= fread(fid3, [2, inf], 'uint8')';
fclose(fid3);
ATRTIME=[];
ANNOT=[];
sa=size(A);
saa=sa(1);
i=1;
while i<=saa
    annoth=bitshift(A(i,2),-2);
    if annoth==59
        ANNOT=[ANNOT;bitshift(A(i+3,2),-2)];
        ATRTIME=[ATRTIME;A(i+2,1)+bitshift(A(i+2,2),8)+...
                bitshift(A(i+1,1),16)+bitshift(A(i+1,2),24)];
        i=i+3;
    elseif annoth==60
        % nothing to do!
    elseif annoth==61
        % nothing to do!
    elseif annoth==62
        % nothing to do!
    elseif annoth==63
        hilfe=bitshift(bitand(A(i,2),3),8)+A(i,1);
        hilfe=hilfe+mod(hilfe,2);
        i=i+hilfe/2;
    else
        ATRTIME=[ATRTIME;bitshift(bitand(A(i,2),3),8)+A(i,1)];
        ANNOT=[ANNOT;bitshift(A(i,2),-2)];
   end;
   i=i+1;
end;
ANNOT(length(ANNOT))=[];       % last line = EOF (=0)
ATRTIME(length(ATRTIME))=[];   % last line = EOF
clear A;
ATRTIME= (cumsum(ATRTIME))/sfreq;
ind= find(ATRTIME <= TIME(end));
ATRTIMED= ATRTIME(ind);
ANNOT=round(ANNOT);
ANNOTD= ANNOT(ind);

%------ DISPLAY DATA ------------------------------------------------------
figure(1); clf, box on, hold on
plot(TIME, M(:,1),'r');
%暂时不要画出通道2的数据
% if nosig==2
%     plot(TIME, M(:,2),'b');
% end;
for k=1:length(ATRTIMED)
    text(ATRTIMED(k),0,num2str(ANNOTD(k)));
end;
xlim([TIME(1), TIME(end)]);
xlabel('Time / s'); ylabel('Voltage / mV');
string=['ECG signal ',DATAFILE];
title(string);
fprintf(1,'\\n$> DISPLAYING DATA FINISHED \n');
%-------------补充对心电的每一周期的平均值计算与处理-----------------------
B=zeros(1,length(ATRTIMED)-1);
for k=1:length(ATRTIMED)-1
    B(k)=ATRTIMED(k+1)-ATRTIMED(k);
end;
Average_internal = mean(B(:));
fprintf('\\n$> ECG周期平均时长为：%8.3f秒 \n',Average_internal);
a1=round(sfreq*Average_internal/2); %一个ECG周期的平均值的采样点数的一半
C=zeros(2*a1+1,length(ATRTIMED)-4); %生成一个2a1行，图中剔除第一个和最后一个峰值的剩余峰值数个列的矩阵
D=zeros(1,length(ATRTIMED)-4); %存放标注的数字，生成一个1行，图中剔除第一个和最后一个峰值的剩余峰值数个列的矩阵
ATRTIMED1=zeros(1,length(ATRTIMED)-4);
for k=1:length(ATRTIMED)-4
    ATRTIMED1(k)=round(ATRTIMED(k+2)*sfreq); %标注处的点是第几个采样点
end;
for k=1:length(ATRTIMED)-4
%     lower_limit = ATRTIMED1(k) - a1
%     upper_limit = ATRTIMED1(k) + a1
    C(:,k)=M(ATRTIMED1(k)-a1:ATRTIMED1(k)+a1,1);
end;

C = C';

for k=1:length(ATRTIMED)-4
    D(1,k)=ANNOTD(k+2);
end;

D = D';
test_normal = C(1:80,:);
test_label_normal = D(1:80,:);
train_normal = C(81:94,:);
train_label_normal = D(81:94,:);
save test_normal test_normal
save test_label_normal test_label_normal
save train_normal train_normal
save train_label_normal train_label_normal

% -------------------------------------------------------------------------
fprintf(1,'\\n$> ALL FINISHED \n');