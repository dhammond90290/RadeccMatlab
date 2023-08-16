function [OutputRow, Header]=RadeccWorkup15(filename)
%version 15 is from 14 (was from RadeccWorkup13) adjusted for max 600 min, has Output (1,41) for quality
%Also adjusts lines where stop and group are read, to hopefully work correctly.
%this program will read text files from Radecc output, calculate chance
% counts for each interval read, subtract these from the interval cpm and
% find corr cpm in each window.
% The corr cpm will be used to calculate a simple average for all
% intervals used.  An option to do weighted average with weighting based on 
% uncertainty in corr cpm has been commented out.  
% Operator can choose range of lines to be included in calculation.
% Range can be different for 219 and 220.
% Data are screened to find anomalous counts based on total count rates >100 cpm
%   and then to find 219 or 220 more than 3.5 ssd from mean (really
%   only works for 219, as ssd for 220 can get large due to increasing cc
%   correction.  Program will ask if each anomalous line should be ignored.
% Program assumes intervals are always equal times for a given file
% THE FILE USED MUST HAVE A FOLDER WITH DATA FILES AND HAVE THE FUNCTION
% LLSQ INCLUDED, AS THIS FUNCTION IS CALLED.

% THIS HAS BEEN ADAPTED TO RUN ON MAC VERSION OF MATLAB.  THERE MAY BE A
% PROBLEM IN RadeccFolder IF IT IS RUN ON PC MATLAB VERSION, DUE TO BACKSLASH CONVENTION.

%THIS VERSION CONTAINS THE NEW ALGORITHMS PROPOSED BY HAMMOND IN JULY 2014

% This version will output results for 220 window after 180 minutes using
% both single line calculation and using cc correction for each line.

%Create a new version in 2022 that uses new Hammond CC corrections.  In
%this version, CC for the 219 and 220 windows will be only dependent on the
%ac counts in their own window. No correction for 2 stage 219 decays in the 
% 220 window are made, as these should occur randomly throughout the count

% settings used for channel cross talk and must be changed if
% gate timing differs from standard.  These are in variables B41, B42, B43
% These depend on half-lives and gate settings
% lambda is decay constant in ms
% t11 and t12 are opening and closing of window 1 in ms
% t21 and t22 are opening and closing of window 2 in ms
% f11 and f12 are fractions of 215Po decaying in windows 1 and 2
% f21 and f22 are fractions of 216Po decaying in windows 1 and 2
t11=0.01;
t12=5.6;
t21=5.61;
t22=600;
lambda215=(log(2)/1.783);
lambda216=(log(2)/145);
f11=exp(-lambda215*t11)-exp(-lambda215*t12);
f12=exp(-lambda215*t21)-exp(-lambda215*t22);
f21=exp(-lambda216*t11)-exp(-lambda216*t12);
f22=exp(-lambda216*t21)-exp(-lambda216*t22);
Det=f11*f22-f12*f21;  
%this is the determinant for the inverse solution to 
% the equations A219=m1*C1-m2*C2 and A220=m3*C2-m4*C1
% where A is detected cpm and C1 and C2 are observed counts in 
%  windows 1 and 2 (corrected for cc and bkg)
m1=f22/(Det);
m2=f21/(Det);
m3=f11/(Det);
m4=f12/(Det);

% We actually want to get detected events in window 1 that belong to 219  
% and events in window two that belong to 220.  Then we need to divide 
% the equations above by either m1 or m3 to find these values
%  This is Giffin correction of approx 2.5% of 220 in 219 channel

B12=m2/m1;
B21=m4/m3;

%CC corrections use constants based on gate timing.  If different
%settings are used, these variables must be changed above at t11 to t22
%tg1 and tg2 below are in minutes.  tg is their sum
    tg1 = (t12-t11)/(60*1000);
    tg2 = (t22-t21)/(60*1000);
    tg = tg1+tg2;

% We also need to know that probability that a 216Po daughter produced by a
% 220Rn parent after the 220 window has opened, will decay during the
% window.  This is found by integrating the function
% (1/tg)(1-exp(-lambda*t)) over time 0 to t22 (in ms).  This parameter is used in
% the correction for non-random coincidence events from the 220 daughter.

pf216 = 1 + (exp(-lambda216*(t22-t21)) - 1)/(lambda216*(t22-t21));

%testcorr is a variable that can be set to determine where to cut off
%counting for 219 or 220.  It is the ratio of chance counts to cpm for
%these windows.  Currently set at 85%, where cc is about 6x corr cpm
testcorr = .85;

%max_minutes is a variable that can be set to ignore lines beyond this time
max_minutes = 900;%was set to 480

%When program screens for noise, it uses ndev*sample std deviation. Set ndev here
ndev=3.3; %this criterion will eliminate lines with <0.3% chance of random occurence

% THIS PROGRAM ASSUMES BACKGROUNDS ARE ZERO, BUT THIS CAN BE CHANGED HERE
BKG219=0; %BKG AT 219 IN cpm
BKG220=0; %BKG AT 220 IN cpm


%First part of code modified from Nick's RadeccRead.m

clear data
clear Header
clear wt

%Added these lines for tests, but do not use
%filename = 'GT_36-16-2A_ch1_20140216c'
%filename = 'Std_IAEA-E_ch1_20140310a'

fid=fopen(filename)

if strcmp(filename(end-3:end),'.txt')==1
    databyline=textread(filename, '%s', 'delimiter', '\n'); %read each line from text fill   
    test = 1
else
    databyline=textscan(fid,'%*c %[^"]'); %command skips (* means skip) first character
    % and reads to double quote "
    databyline=databyline{1,1}; %break nested cells into just cells with lines of data (every other row is empty)
    test = 2
end

%remove empty rows
for i=1:length(databyline(:,1))
    ltemp(i,1)=length(databyline{i,1});
end
databyline=databyline(ltemp>=3,1);





lines=0;  %this variable will be number of data lines read in
        %the array data has values read in from each line
        % the array wt will be weighting for line, first set to 1.
shortform=0;
for i=2:length(databyline) %skip first row as it could by chance have 46 or 49 chars
    thisline=textscan(databyline{i,1},'%c');
    if length(thisline{1,1})==49 ||  length(thisline{1,1})==46 %read 46 or 49 length rows as data rows
        lines=lines+1;
        j=lines;
        a=textscan(databyline{i,1},'%f');
        data(j,:)=a{1,1}'; wt(j) = 1;
        if length(thisline{1,1})==46
            shortform=1; %if length of data line is 46, flag as shortform to read stop/start differntly
        end
    end
    
end
disp(' Lines read =  ')
lines=lines-1 %DH added this 7/29/17
fclose(fid);


% Make sure file has at least 3 lines and if not, return and get next file
if lines<=2
    disp('This file had 2 or less lines... skipping that garbage')
  % the return statement exits the function and goes back to next data file
    OutputRow{1,1}=filename; %filename
    temp=textscan(databyline{2,1},'%10c %s %s'); %start date time strings
    OutputRow{1,2}=temp{1,2}{1,1}; %start date
    OutputRow{1,3}=temp{1,3}{1,1}; %start time
    
    temp=textscan(databyline{lines+4,:},'%15c %s %s'); %stop date time stings
    OutputRow{1,4}=temp{1,2}{1,1}; %stop date
    OutputRow{1,5}=temp{1,3}{1,1}; %stop time
    return
end

% data is a file with the lines read from the Radecc file
% It has columns of time, 219cpm,219counts for interval, 220cpm, 220counts
% per interval, totalccpm, total counts per interval. Time and cpm are integrated values since count start
% Number of lines in data = number of blocks read out



%Now zero del_t, 219, 220, total All are 1-D arrays, length will be = data length
% set these arrays at zero to start
% xcalc array will be filled with data to plot: line#,obs219cpm,sig219,obs220cpm,sig220,tot cpm
% These obs219cpm and obs220cpm will be corrected for chance counts, data for each line

del_t = [ ]; %zeros(500);
del_219 = [ ] ;%zeros(500);
del_220 = [ ]; %zeros(500);
del_tot = [ ] ;%zeros(500);
xcalc = [ ] ; %zeros(500,6);



%del_t is an array with interval length, del_219 etc is number of counts for interval i
%lines = j and  j is the number of lines of data in the array

del_t(1) = data(1,1)
del_219(1) = data(1,3);
del_220(1)= data(1,5);
del_tot(1)= data(1,7);


%if file read for longer than max_minutes, only use counts up to that point
maxlines=round(max_minutes/del_t(1));
if lines>maxlines
    data=data(1:maxlines,:);
    lines=maxlines;
end


for i=2:lines
    del_t(i) = data(i,1)-data(i-1,1);
    del_219(i)= data(i,3);
    del_220(i)= data(i,5);
    del_tot(i)= data(i,7);
end



%FOR EACH LINE, COMPUTE CPM, CHANCE COUNTS, COUNTS CORRECTED FOR CHANCE, ratio of CC to raw cpm
CORR219 = [];
CORR220 = [];
sig219 = [];
sig220 = [];
CPM219 = [];
CPM220 = [];
CPMtot = [];
CC219 = [];
CC220 = [];
calccq219= [];
calccpm219=[];
frac219=[];
calccq220= [];
calccpm220=[];
frac220=[];


for i=1:lines
    
    CPM219(i)=(del_219(i))/del_t(i); 
    CPM220(i)=(del_220(i))/del_t(i);
    CPMtot(i)=(del_tot(i))/del_t(i);
    
    A219=CPMtot(i)-CPM219(i);
    A220=CPMtot(i)-CPM220(i);
    if(A219==0);%Test for zero and set to low value to avoid dividing by zero
        A219=0.0001;
    end
    if(A220==0);%Test for zero and set to low value to avoid dividing by zero
        A220=0.0001;
    end
   
    %2022 program has adjusted these CC calcs below
    CC219(i)=(A219^2)*tg1/(1- A219*tg1);
    CC220(i)=(A220^2)*tg2/(1- A220*tg2);
    %CC220(i)=CC220(i)+((1/f11)*(CPM219(i)-CC219(i)) + (pf216/f22)*(CPM220(i)-CC220(i)))*CC220(i)/A;
    
    % find corrected cpm and sig in cpm for each window.  
    % channel crosstalk will be corrected later
    % uncertainty for each interval is computed but would only be used for
    % version with weighted average
    
    CORR219(i)=(CPM219(i)-BKG219-CC219(i));
    CORR220(i)=(CPM220(i)-BKG220-CC220(i));
    sig219(i) = sqrt(CPM219(i) + CC219(i) + BKG219)/sqrt(del_t(i));
    sig220(i) = sqrt(CPM220(i) + CC220(i) + BKG220)/sqrt(del_t(i));
    
    % xcalc Array below will be used for plots in output
    % xcalc array will be filled with data to plot: line#,obs219cpm,sig219,obs220cpm,sig220,tot cpm
    % These 219cpm and 220cpm have been corrected for chance counts, data for each line
   xcalc(i,1) = i;
    xcalc(i,2) = CORR219(i);
    xcalc(i,3) = sig219(i);
    xcalc(i,4) = CORR220(i);
    xcalc(i,5) = sig220(i);
    xcalc(i,6) = CPMtot(i);
    
    
end

% and display results with sample std deviations
disp('Preliminary aves')
ave219 = mean(CORR219);
stdev219 = std(CORR219);
ave220 = mean((CORR220));
stdev220 = std(CORR220);

%Find line where CC becomes 65% of raw cpm (test criterion can be changed by setting
%variable testcorr to a different value (statement is near beginning of program
%To do this, fit ratio of cc/rawcpm vs line# to a quadratic and get root, where cc/cpm are about to exceed testcorr
%If you want more lines, this can be set after data is plotted.
%Coeffs for quadratic from polyfit should be y=a1*x^2 + a2*x + a3

fitcq219=polyfit(xcalc(:,1),CC219(:),2);
fitcpm219=polyfit(xcalc(:,1),CPM219(:),2);
calccq219=(xcalc(:,1).^2*fitcq219(1)+xcalc(:,1).*fitcq219(2)+fitcq219(3));
calccpm219=(xcalc(:,1).^2*fitcpm219(1)+xcalc(:,1).*fitcpm219(2)+fitcpm219(3));
frac219=calccq219./calccpm219;
% frac219 is a variable with the ratio for each interval of chance counts/cpm for 219 based on fitted curves
% Cannot do this for raw counts as there may be zero cpm for some intervals

fitcq220=polyfit(xcalc(:,1),CC220(:),2);
fitcpm220=polyfit(xcalc(:,1),CPM220(:),2);
calccq220=(xcalc(:,1).^2*fitcq220(1)+xcalc(:,1).*fitcq220(2)+fitcq220(3));
calccpm220=(xcalc(:,1).^2*fitcpm220(1)+xcalc(:,1).*fitcpm220(2)+fitcpm220(3));
frac220=calccq220./calccpm220;
% frac220 is like frac219

% fit linear regression to total cpm
fitcpmtot=polyfit(xcalc(:,1),CPMtot(:),1);

% find where fit to line has (total cpm - corr219 - corr220 cpm) = 15 cpm
lines_totcpm15=round((15+ave219+ave220-fitcpmtot(2))/fitcpmtot(1));
if lines_totcpm15<3
    lines_totcpm15=lines;
end

%Now Test to see if ratio exceeds test value
for i=1:lines
    if frac219(i)>testcorr
        lines_test219=i; 
        break
    else
        lines_test219=lines;
    end
end
    
%Must Check to be sure enough lines (3 or more) are included
%If not, all lines are included and operator will hand pick interval to use
if lines_test219<3
    disp('cc/cpm fraction of 219 array is ')
    disp(frac219)
    disp('Not enough good 219 lines. Set 219 lines to max lines')
    lines_test219=lines
end
    
    for i=1:lines
    if frac220(i)>testcorr
        lines_test220=i; 
        break
    else
        lines_test220=i;
    end
end
        
if lines_test220<3
    disp('cc/cpm = fraction for 220 array is ')
    disp(frac220)
    disp('Good 220 lines<3. Set 220 lines where total cpm less219&220=15')
    lines_test220=lines_totcpm15
    if lines_totcpm15>lines
        lines_test220=lines
    end
end


%Display averages, but not sample standard deviations
disp(['ave219 = ', num2str(mean(CORR219))])
%disp(['stdev219 = ', num2str(std(CORR219))])
disp(['ave220 = ', num2str(mean((CORR220)))])
%disp(['stdev220 = ', num2str(std(CORR220))])

% Plot data vs line number for selected interval.  Use max lines (set at up to 72) for total. 
%Could display equations for lines here by removing comments below
x1 = 1;
x2 = lines;
x1_tot =x1; x2_tot =x2; 
x1_219 =x1; x2_219 =lines_test219;
x1_220 =x1; x2_220 =lines_test220;
        close all
        figure('name',filename)
        subplot(3,1,1)
        plot((x1_tot:x2_tot),CPMtot(x1_tot:x2_tot));
        title('Corrected CPM vs Line #');
        hold on
        coefftot=regstats(CPMtot(x1_tot:x2_tot),(x1_tot:x2_tot),'linear',{'beta' 'rsquare'});
        plot((x1_tot:x2_tot),(x1_tot:x2_tot)*coefftot.beta(2)+coefftot.beta(1),'g')
%        disp(['total fit:  y = ',num2str(coefftot.beta(2)),'x + ',num2str(coefftot.beta(1))])
%        disp(['rsquared = ',num2str(coefftot.rsquare(1))])
        subplot(3,1,2)
        plot((x1_219:x2_219),CORR219(x1_219:x2_219));
        title('Corrected 219 CPM vs Line #');
        hold on
        coeff219=regstats(CORR219(x1_219:x2_219),(x1_219:x2_219),'linear',{'beta' 'rsquare'});
        plot((x1_219:x2_219),(x1_219:x2_219)*coeff219.beta(2)+coeff219.beta(1),'g')
%        disp(['219 fit:  y = ',num2str(coeff219.beta(2)),'x + ',num2str(coeff219.beta(1))])
%        disp(['rsquared = ',num2str(coeff219.rsquare(1))])
        subplot(3,1,3)
        plot((x1_220:x2_220),CORR220(x1_220:x2_220));
        title('corrected 220cpm vs Line#')
        xlabel('Line #');
        hold on
        coeff220=regstats(CORR220(x1_220:x2_220),(x1_220:x2_220),'linear',{'beta' 'rsquare'});
        plot((x1_220:x2_220),(x1_220:x2_220)*coeff220.beta(2)+coeff220.beta(1),'g')
%        disp(['220 fit:  y = ',num2str(coeff220.beta(2)),'x + ',num2str(coeff220.beta(1))])
%        disp(['rsquared = ',num2str(coeff220.rsquare(1))])

%Now select line ranges to use for calculation and to replot data.
%User can choose interval to use in calculation of aves and slopes or use
%      default values already found
%If a value out of range is chosen, it assigns value of lines
for cycles=1:20 
    disp(' ')
    disp('Enter 0 to change all CPM ranges')
    disp('Enter 1 to change total CPM range')
    disp('Enter 2 to change 219 CPM range')
    disp('Enter 3 to change 220 CPm range')
    flag=input('or enter 5 to continue & plot final output  :');
    
    if flag==5
        break
    end
    
    if flag==1;
        x1_tot = input('enter new x1 value for tot cpm ');
        x2_tot = input('enter new x2 value for tot cpm ');
        if x2_tot>lines || x1_tot<1 || x1_tot>x2_tot
            x2_tot=x2;
            x1_tot=x1;
        end
 
    end
    
    if flag==2
        x1_219 = input('enter new x1 value for 219  ');
        x2_219 = input('enter new x2 value for 219 ');
        if x2_219>lines || x1_219<1 || x1_219>x2_219
            x2_219=x2;
            x1_219=x1;
        end
    end
    
    if flag==3
        x1_220 = input('enter new x1 value for 220  ');
        x2_220 = input('enter new x2 value for 220  ');
        if x2_220>lines || x1_220<1 || x1_220>x2_220
            x2_220=x2;
            x1_220=x1;
        end
    end
    
    if flag==0
        x1_220 = input('enter new x1 value for all  ');
        x2_220 = input('enter new x2 value for all  ');
        if x2_220>lines || x1_220<1 || x1_220>x2_220
            x2_220=x2;
            x1_220=x1;
        end
        x2_219=x2_220;
        x2_tot=x2_220;
        x1_219=x1_220;
        x1_tot=x1_220;
    end
    
    close all

%Plot Data    
    figure('name',filename)
    subplot(3,1,1)
    plot((x1_tot:x2_tot),CPMtot(x1_tot:x2_tot));
    title('Corrected CPM vs Line #');
    hold on
    coefftot=regstats(CPMtot(x1_tot:x2_tot),(x1_tot:x2_tot),'linear',{'beta' 'rsquare'});
    plot((x1_tot:x2_tot),(x1_tot:x2_tot)*coefftot.beta(2)+coefftot.beta(1),'g')
    disp(['total fit:  y = ',num2str(coefftot.beta(2)),'x + ',num2str(coefftot.beta(1))])
    disp(['rsquared = ',num2str(coefftot.rsquare(1))])
    subplot(3,1,2)
    plot((x1_219:x2_219),CORR219(x1_219:x2_219));
    title('Corrected 219 CPM vs Line #');
    hold on
    coeff219=regstats(CORR219(x1_219:x2_219),(x1_219:x2_219),'linear',{'beta' 'rsquare'});
    plot((x1_219:x2_219),(x1_219:x2_219)*coeff219.beta(2)+coeff219.beta(1),'g')
    disp(['219 fit:  y = ',num2str(coeff219.beta(2)),'x + ',num2str(coeff219.beta(1))])
    disp(['rsquared = ',num2str(coeff219.rsquare(1))])
    subplot(3,1,3)
    plot((x1_220:x2_220),CORR220(x1_220:x2_220));
    title('corrected 220cpm vs Line#')
    xlabel('Line #');
    hold on
    coeff220=regstats(CORR220(x1_220:x2_220),(x1_220:x2_220),'linear',{'beta' 'rsquare'});
    plot((x1_220:x2_220),(x1_220:x2_220)*coeff220.beta(2)+coeff220.beta(1),'g')
    disp(['220 fit:  y = ',num2str(coeff220.beta(2)),'x + ',num2str(coeff220.beta(1))])
    disp(['rsquared = ',num2str(coeff220.rsquare(1))])
    
end

% Scan selected data range to use to see if any line has anomalous counts for 219 or 220. 

%FIRST check for noise by looking for total count rate intervals >100 cpm

noisylines = 0;

disp(' ')
disp('....Looking for bad data')
wt=zeros(length(CORR219),2)+1;
for i=min(x1_219,x1_220):max(x2_219,x2_220)
    if CPMtot(i)>80
       disp(['line number ', num2str(i), ' looks bad as cpmtot>80']) 
       value=CPMtot(i)
       wt(i,1) = input( ' enter 0 to ignore line or 1 to keep line  :');
       wt(i,2) = wt(i,1);
       noisylines = noisylines + (1-wt(i,1));
    end
end


%  Second average lines and display results with sample std deviations

sumwt219=0;
sum219=0;
sumsq219=0;

for i=x1_219:x2_219
    sum219=sum219+wt(i,1)*CORR219(i);
    sumsq219=sumsq219+wt(i,1)*CORR219(i)*CORR219(i);
    sumwt219=sumwt219+wt(i,1);
end

ave219 = sum219/sumwt219;
sstdev219=sqrt((sumsq219-sumwt219*ave219*ave219)/(sumwt219-1));

sumwt220=0;
sum220=0;
sumsq220=0;

for i=x1_220:x2_220
    sum220=sum220+wt(i,2)*CORR220(i);
    sumsq220=sumsq220+wt(i,2)*CORR220(i)*CORR220(i);
    sumwt220=sumwt220+wt(i,2);
 end

ave220 = sum220/sumwt220;
sstdev220=sqrt((sumsq220-sumwt220*ave220*ave220)/(sumwt220-1));

disp(' ')
disp('Using the selected sections minus really bad noise:')
disp('NOTE - these are sample standard deviations')
disp(['ave219 = ', num2str(ave219)])
disp(['sstdev219 = ', num2str(sstdev219)])
disp(['ave220 = ', num2str(ave220)])
disp(['sstdev220 = ', num2str(sstdev220)])


% % % UNUSED CODE BELOW AS ALTERNATE WAY TO DO AVE
% % %make matrix to NaN out bad values found earlier
% % % temp= wt==0;
% % % nanwt=wt;
% % % nanwt(temp)=NaN;
% % % nanwt=nanwt(:,1);
% % % 
% % % ave219wt = nanmean(CORR219.*nanwt');
% % % stdev219wt = nanstd(CORR219.*nanwt');
% % % ave220wt = nanmean((CORR220.*nanwt'));
% % % stdev220wt = nanstd(CORR220.*nanwt');


%look to see if counts for any line are more than 3.5 sig from mean, if the
%    interval has more than one count. The selection criterion
%can be modified to change deviations picked by changing variable ndev.
% Program will set weight for these lines to zero, unless over-ridden below
% wt(i,j) array has weighting for 219 (j=1) or 220 (j=2) for each line
% Set weighting to 1 (default) or 0 if counts are anomalous in either

%NOW LOOK FOR OUTLIERS IN 219 AND 220.  CRITERION WAS SET IN BEGINNING, VARIABLE ndev
disp(['...now checking counts outside of ', num2str(ndev),' sig from mean'])
for i=x1_219:x2_219
    test219=(CORR219(i)-ave219)/(ndev*sstdev219);
    if test219>=1 && CORR219(i)*del_t(i)>1 && wt(i,1)==1
        disp(['line number ', num2str(i), ' looks bad for 219: #ssdev is '])
        value=(CORR219(i)-ave219)/sstdev219
        wt(i,1) = input( ' enter 0 to ignore line or 1 to keep line  :');
        wt(i,2) = wt(i,1);
       noisylines = noisylines + (1-wt(i,1));
    end
end
for i=x1_220:x2_220
    test220=(CORR220(i)-ave220)/(ndev*sstdev220);
    if test220>=1 && CORR220(i)*del_t(i)>1 && wt(i,1)==1
        disp(['line number ', num2str(i), ' looks bad for 220: #ssdev is '])
        value=(CORR220(i)-ave220)/sstdev220
        wt(i,1) = input( ' enter 0 to ignore line or 1 to keep line  :');
        wt(i,1) = wt(i,2);
        noisylines = noisylines + (1-wt(i,1));
   end
end
disp('...check done')

%re-make matrix to NaN out bad values found earlier
temp= wt==0;
nanwt=wt;
nanwt(temp)=NaN;
nanwt=nanwt(:,1);

%Check again to make sure there are at least 3 good lines left
if lines<=2
    disp('This file had 2 or less lines... skipping that garbage')
  % the return statement exits the function and goes back to next data file
    OutputRow{1,1}=filename; %filename
    temp=textscan(databyline{2,1},'%10c %s %s'); %start date time strings
    OutputRow{1,2}=temp{1,2}{1,1}; %start date
    OutputRow{1,3}=temp{1,3}{1,1}; %start time
    
    temp=textscan(databyline{lines+4,:},'%15c %s %s'); %stop date time stings
    OutputRow{1,4}=temp{1,2}{1,1}; %stop date
    OutputRow{1,5}=temp{1,3}{1,1}; %stop time
    return
end



%Do linear least squares fit using function LLSQ. Itreturns 2x2 array with line 1 of m and sig m, line 2 of b and sig b

disp('New fit ignoring bad data')
mbtot=LLSQ([(1:lines)',CPMtot',wt(:,1)],x1_tot,x2_tot);
disp(['Final Tot fit:  y = (' , num2str(mbtot(1,1)) ,'+/-', num2str(mbtot(1,2)) ,')x  +  (' , num2str(mbtot(2,1)) , '+/-', num2str(mbtot(2,2)), ')' ])
mb219=LLSQ([(1:lines)',CORR219',wt(:,1)],x1_219,x2_219);
disp(['Final 219 fit:  y = (' , num2str(mb219(1,1)) ,'+/-', num2str(mb219(1,2)) ,')x  +  (' , num2str(mb219(2,1)) , '+/-', num2str(mb219(2,2)), ')' ])
mb220=LLSQ([(1:lines)',CORR220',wt(:,2)],x1_220,x2_220);
disp(['Final 220 fit:  y = (' , num2str(mb220(1,1)) ,'+/-', num2str(mb220(1,2)) ,')x  +  (' , num2str(mb220(2,1)) , '+/-', num2str(mb220(2,2)), ')' ])




% % % % Find weighted means for 219 and 220 windows FOR VALUES SET AS x1 and x2
% % % % Possible problem is whether at low count rates, wt is too heavy in
% % % % non-zero channels.  THIS OPTION IS IGNORED FOR NOW.
% % % sum219 = 0.0;
% % % sum220 = 0.0;
% % % sumwt219 = 0.0;
% % % sumwt220 = 0.0;
% % % nskip = 0;
% % % for i=x1:x2
% % %
% % %     wt219=1/(sig219(i))^2;
% % %     if wt(i,1)<0.1
% % %         wt219=0;
% % %         nskip=nskip+1;
% % %     end
% % %
% % %     wt220=1/(sig220(i))^2;
% % %     if wt(i,2)<0.1
% % %         wt220=0;
% % %     end
% % %
% % %     sum219 = sum219 + wt219*CORR219(i);
% % %     sum220 = sum220 + wt220*CORR220(i);
% % %     sumwt219 = sumwt219 + wt219;
% % %     sumwt220 = sumwt220 + wt220;
% % % end
% % % AveCPM219 = sum219/sumwt219
% % % AveCPM220 = sum220/sumwt220
% % % sigAve219 = sqrt(1/sumwt219)
% % % sigAve220 = sqrt(1/sumwt220)



%NOW compute simple averages for x1 to x2
ssum219=0.0;
ssum220 =0.0;
ssumwt219 =0.0;
ssumwt220 =0.0;
ssq219 = 0.0;
ssq220 = 0.0;
for i=x1_219:x2_219
    wt219= wt(i,1);
    ssumwt219 = ssumwt219 + wt219;
    ssum219 = ssum219 + wt219*CORR219(i);
    ssq219=ssq219 + wt219*CORR219(i)^2;
end
for i=x1_220:x2_220
    wt220= wt(i,2);
    ssumwt220 = ssumwt220 + wt220;
    ssum220 = ssum220 + wt220*CORR220(i);
    ssq220=ssq220 + wt220*CORR220(i)^2;
end
sAveCPM219 = ssum219/ssumwt219; %simple avg
sAveCPM220 = ssum220/ssumwt220;  %simple avg
ssigAve219 = sqrt((ssq219-ssumwt219*(sAveCPM219)^2)/(ssumwt219*(ssumwt219-1)));  %sdom for simple avg
ssigAve220 = sqrt((ssq220-ssumwt220*(sAveCPM220)^2)/(ssumwt220*(ssumwt220-1)));  %sdom for simple avg 

disp(' ')
disp(['Simple Avg 219 window cpm = ', num2str(sAveCPM219)])
disp(['Simple Avg 219 window sigma = ', num2str(ssigAve219)])
disp(['Simple Avg 220 window cpm = ', num2str(sAveCPM220)])
disp(['Simple Avg 220 window sigma= ', num2str(ssigAve220)])

% Now calc data for 220 after 180 minutes with cc correction for 180m ave
% (linear model) or using cc for each line (nonlinear)

% First do calculation of counts after 180 min without cc corrections
% lines180 is variable with number of lines to get to 180 minutes

lines180=180/del_t(1) + 0.2; % this should make lines run through 180 min if first line is a little too big.

ssum219 = 0.0;
ssum220 =0.0;
ssumtime =0.0;
ssumtot = 0.0;

if x2_220>lines180  %check number of lines that are good up to 180 min
    lines220=lines180;
else
    lines220=x2_220;
end
for i=x1_220:lines220
    wt220= wt(i,2);
    ssumtime = ssumtime + wt220*del_t(i);
    ssum219 = ssum219 + wt220*del_219(i);
    ssum220 = ssum220 + wt220*del_220(i);
    ssumtot = ssumtot + wt220*del_tot(i);
end

% Second, calc chance counts for this group with ave cpms for full 180 min
% Variables tg1 and tg2 were set at beginning of program, depending on gate timing
    linCPM219=(ssum219)/ssumtime; 
    linCPM220=(ssum220)/ssumtime;
    linCPMtot=(ssumtot)/ssumtime;
    
    A=linCPMtot-linCPM219-linCPM220;
    
    linCC219=(tg1/tg)*(A^2)*tg/(1- A*tg);
    linCC220=(tg2/tg)*(A^2)*tg/(1- A*tg);
    linCC220=linCC220+((1/f11)*(linCPM219-linCC219) +(pf216/f22)*(linCPM220-linCC220))*linCC220/A;

     
    % find corrected cpm and sig in cpm for each window.  
    % channel crosstalk will be corrected later
    % uncertainty for each interval is computed but would only be used for
    % version with weighted average
    
    linCORR219=(linCPM219-BKG219-linCC219);
    linCORR220=(linCPM220-BKG220-linCC220);
    linsig219 = sqrt(linCPM219 + linCC219 + BKG219)/sqrt(ssumtime);
    linsig220 = sqrt(linCPM220 + linCC220 + BKG220)/sqrt(ssumtime);


% Third, use up to 180 minutes with line by line chance correction (non-linear calc) for 220 calc

ssum220 =0.0;
ssumwt220 =0.0;
ssq220 = 0.0;

for i=x1_220:lines220
    wt220= wt(i,2);
    ssumwt220 = ssumwt220 + wt220;
    ssum220 = ssum220 + wt220*CORR220(i);
    ssq220=ssq220 + wt220*CORR220(i)^2;
end
ave220cpm180nonlin = ssum220/ssumwt220; %nonlinear avg sigma at 180 min
sig220Ave180nonlin = sqrt((ssq220-ssumwt220*(sAveCPM220)^2)/(ssumwt220*(ssumwt220-1)));  



% Next few statements correct counts in each window for cross talk of
% windows.  Variables needed to get meas cpm in each window are B12,B21
% These Constants depend on gate settings and isotope decay
% constants and were set at beginning of program
% Calculated cpm for stds and samples is observed cpm, 
%   not corrected for the 9-12% of decays not captured by windows
% These are the final numbers desired and sigma is computed based on sdom 
%   for all intervals used

finalCORR219=( sAveCPM219 - B12*sAveCPM220 );
finalCORR220=( sAveCPM220 - B21*sAveCPM219 );
finalsig219 = sqrt(ssigAve219^2 + (B12*ssigAve220)^2);
finalsig220 = sqrt(ssigAve220^2 + (B21*ssigAve219)^2);

final220cpm180nonlin=( ave220cpm180nonlin - B21*sAveCPM219 );
final220sig180nonlin = sqrt(sig220Ave180nonlin^2 + (B21*ssigAve219)^2);

final219cpm180lin=( linCORR219 - B12*linCORR220 );
final220cpm180lin=( linCORR220 - B21*linCORR219 );
final219sig180lin = sqrt(linsig219^2 + (B12*linsig220)^2);
final220sig180lin = sqrt(linsig220^2 + (B21*linsig219)^2);


disp(' ')
disp(['Final value 219 cpm = ', num2str(finalCORR219)])
disp(['Final value 219 cpm sig = ', num2str(finalsig219)])
disp(['Final value 220 cpm = ', num2str(finalCORR220)])
disp(['Final value 220 cpm sig= ', num2str(finalsig220)])
disp(' ')

%Now compute the net counts used in the computation
net219cnts = ssumwt219*del_t(1)*sAveCPM219;
net220cnts180 = ssumwt220*del_t(1)*ave220cpm180nonlin;


%NOW PLOT FINAL DATA USED
close all
    figure('name',filename)
    subplot(3,1,1)
    plot((x1_tot:x2_tot),nanwt(x1_tot:x2_tot).*CPMtot(x1_tot:x2_tot)');
    title('FINAL Corrected CPM vs Line #');
    hold on
    plot((x1_tot:x2_tot),(x1_tot:x2_tot)*mbtot(1,1)+mbtot(2,1),'r')

    subplot(3,1,2)
    plot((x1_219:x2_219),nanwt(x1_219:x2_219).*CORR219(x1_219:x2_219)');
    title('FINAL Corrected 219 CPM vs Line #');
    hold on
    plot((x1_219:x2_219),(x1_219:x2_219)*mb219(1,1)+mb219(2,1),'r')

    subplot(3,1,3)
    plot((x1_220:x2_220),nanwt(x1_220:x2_220).*CORR220(x1_220:x2_220)');
    title('FINAL corrected 220 CPM vs Line #')
    xlabel('Line #');
    hold on
    plot((x1_220:x2_220),(x1_220:x2_220)*mb220(1,1)+mb220(2,1),'r')





Header={'Filename', 'DateStarted', 'TimeStarted', 'Stop Date', 'Stop Time' 'Channel#', ...
    'Total Intervals', 'Noise lines', 'Lines 219','Lines 220',...
    'Final Corrected 219', 'fin cor 219 sig', 'Final Corrected 220', 'fin cor 220 sig',  ...
    'slope219 %/hr', '219slope sig/slope','int219 cpm','int 219 error','Simple Ave 219cpm','simp avg 219 sig', 'slope220 %/hour', '220sig slope/slope', 'int220cpm', 'int 220 error', ...
    'Simple Average 220','simp avg 220 sig', 'totcpmslope', 'slope totcpm error','totcpmint', 'int totcpm error', 'net 219 counts', 'net 220 counts at180m nonlin', ...
    'Matlab Datenum StartTime', 'Matlab Datenum EndTime', '180m linfinal219cpm', '180m linfinal219cpmsig', ...
    '180m linfinal220cpm', '180m linfinal220cpmsig', '180m nonlinfinal220cpm', '180m nonlinfinal220cpmsig', 'Quality'};

%Make output row
%For slopes of 219 and 220, output is percent change per hour
%For slope uncertainties of 219 and 220, output is ratio of uncertainty to slope
%Line slope and intercept coeffs are not corrected for cross talk, so normalize with simple ave for 219 and 220
%For total cpm, slope and intercepts are in cpm per line

%find number of lines per hour
lines60=60/del_t(1);

OutputRow{1,1}=filename; %filename

OutputRow{1,1}=filename; %filename

if shortform==0 %if line length was normal read dates this way
    temp=textscan(databyline{2,1},'%10c %s %s'); %start date time strings
    OutputRow{1,2}=temp{1,2}{1,1}; %start date
    OutputRow{1,3}=temp{1,3}{1,1}; %start time
    
    temp=textscan(databyline{lines+4,:},'%15c %s %s'); %stop date time stings
    OutputRow{1,4}=temp{1,2}{1,1}; %stop date
    OutputRow{1,5}=temp{1,3}{1,1}; %stop time
end

if shortform==1 %if line length was short read dates this way
    temp=textscan(databyline{2,1},'%10c %s %s'); %start date time strings
    OutputRow{1,2}=temp{1,2}{1,1}; %start date
    OutputRow{1,3}=temp{1,3}{1,1}; %start time
    
    temp=textscan(databyline{lines+4,:},'%15c %s %s'); %stop date time stings
    %DH changed from lines +3 to lines+4
    OutputRow{1,4}=temp{1,2}{1,1}; %stop date
    OutputRow{1,5}=temp{1,3}{1,1}; %stop time
end
    
%temp=textscan(filename, '%[^_] %*1c %[^_] %*1c %[^_] %*1c %[^_]  %*1c %[^_]'); %read filename in parts seperated by _
%chan=textscan(databyline{end-2,:},'%*s %*s %f');OLD, REPLACED BY FOLLOWING
chan=textscan(databyline{lines+9,:},'%*s %*s %f');% DH changed from 9 to 9
OutputRow{1,6}=chan{1,1}; %channel #

OutputRow{1,7}=lines; %total intervals
OutputRow{1,8}=noisylines; %intervals ignored for noise
OutputRow{1,9}=ssumwt219; %intervals used for 219
OutputRow{1,10}=x2_220 - noisylines; %intervals used for 220

OutputRow{1,11}=finalCORR219; %final corrected averages and uncertainties
OutputRow{1,12}=finalsig219;
OutputRow{1,13}=finalCORR220;
OutputRow{1,14}=finalsig220;

OutputRow{1,15}=100*lines60*(mb219(1,1)/sAveCPM219); %219 slope in % per hour
OutputRow{1,16}=mb219(1,2)/mb219(1,1); %ratio of slope err/slope for 219
OutputRow{1,17}=mb219(2,1); %219 int
OutputRow{1,18}=mb219(2,2); %219 int err
OutputRow{1,19}=sAveCPM219; %219 simple average
OutputRow{1,20}=ssigAve219; %219 simple average err

OutputRow{1,21}=100*lines60*(mb220(1,1)/sAveCPM220); %220 slope in % per hour
OutputRow{1,22}=mb220(1,2)/mb220(1,1); %ratio of slope err/slope for 220
OutputRow{1,23}=mb220(2,1); %220 int
OutputRow{1,24}=mb220(2,2); %220 int err
OutputRow{1,25}=sAveCPM220; %220 simple average
OutputRow{1,26}=ssigAve220; %220 simple average err

OutputRow{1,27}=mbtot(1,1); %Tot slope
OutputRow{1,28}=mbtot(1,2); %Tot slope err
OutputRow{1,29}=mbtot(2,1); %Tot int
OutputRow{1,30}=mbtot(2,2); %Tot int err

OutputRow{1,31}=net219cnts;
OutputRow{1,32}=net220cnts180;

OutputRow{1,33}=datenum(strcat(OutputRow(1,2),'_',OutputRow(1,3))); %days since beginning of 01/01/0000 and start
OutputRow{1,34}=datenum(strcat(OutputRow(1,4),'_',OutputRow(1,5))); %days since beginning of 01/01/0000 and end

OutputRow{1,35}=final219cpm180lin;
OutputRow{1,36}=final219sig180lin;
OutputRow{1,37}=final220cpm180lin;
OutputRow{1,38}=final220sig180lin;

OutputRow{1,39}=final220cpm180nonlin; %corrected 220 after 180 min (multi-line)
OutputRow{1,40}=final220sig180nonlin; %sig in corrected 220 after 180 min
    disp('Enter quality')
    quality = input(' 1=good, 2=noise, 3=drop>10% over range used, 4=increase>20% over range ');
OutputRow{1,41} = quality;
 
%FOREXCEL=[Header;OutputRow];







