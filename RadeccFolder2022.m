%function [Output]=RadeccFolder(foldername)
%this function takes a folder in the current matlab directory containing
%Radecc output files, runs RadeecWorkupXX routine on each file and
%creates an excel file with the results
%You must have RadeccWorkup##, LLSQ, and the data folder all in one folder
%You will want to periodically copy results from the variable Output and
%paste them to an excel file
%You need to first create a file 'RadeccCalcArchive' to store results
%ON LINE 35, SET THE VERSION OF RadeccWorkup THAT YOU WANT TO USE.

%The next line will identify the file RadeccCalcArchive where data is filed
%THIS OPTION IS NOTWORKING 
%      savefile = 'RadeccCalcArchive'

%If running on a PC, pay attention to comment on line 27
nrecord=1;
disp('Enter the foldername in single quotes')
foldername=input('    :')
disp('Enter your starting DATA FILE number')
startnumber=input('(3=first file in folder)  :');
disp('Enter the desired starting row number to write Output in')
% This allows you to use more than one line for multiple work-ups of same
% data file or to continue adding to an existing set of outputs
% Remember that if the program crashes, you lose all calcs not saved to excel
OutputLineNum=input('(2=first row below Header)  :');
directory=dir(foldername); %creates an array with the file names/info



filenumber=startnumber; % this variable keeps track of line# in data file

while filenumber<=length(directory)

    filename=[foldername, '/',directory(filenumber).name]; %may need to change direction of this slash for mac to PC
    [Output(OutputLineNum,:),Header]=RadeccWorkup2022new(filename);
    %Output(filenumber-1,1)={directory(filenumber).name;};
    
    %The next line should save the calculation in a file in case of crash.
    %THIS OPTION IS NOTWORKING. but Nick has new command that overwrites.
    % save(savefile,'Output(OutputLineNum,:)','-ASCII','-append')
     save([foldername,'.mat'],'Output','-mat')
  
    
    Output(1,:)=Header;

    disp('Enter 1-8 to move onto next file')
    flag=input('Enter 9 to re-do same file and add another Output Row   :');
    if flag==9
        filenumber=filenumber-1;
    end
    
    filenumber=filenumber+1;
    OutputLineNum=OutputLineNum+1;
    
    nrecord=nrecord+1;
    if nrecord==15
        disp('you have done 15 runs.  Stop and save variable Output')
        nrecord=1
    end
end


disp('The Matlab variable Output has the data, Copy and Paste into Excel')
