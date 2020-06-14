%Kappa distribution analysis

%% Predefine Paths and Patient ID
% clear variables
clear
clc

%load list of selected patient number
patientnumber = readtable('patientnumber.csv');
patientnumber = table2array(patientnumber);
patientnumber_XsensN = patientnumber(:,1);
patientnumber_XsensC = patientnumber(:,2);
patientnumber_Apple = patientnumber(:,3);

%Define input and output pathways
pathEY ='D:\Box Sync\Box Sync\Elton Shared\All task classifications (E.Y.) V3\';
pathMA ='D:\Box Sync\Box Sync\Elton Shared\All task classifications (M.A.) V3\';
pathTW ='D:\Box Sync\Box Sync\Elton Shared\Natural Xsens classifications (T.W.) V3 csvformat\';
pathOUT = 'D:\Box Sync\Box Sync\Elton Shared\KappaFigures\Summated CM\';

%load list of files in each rater's folder for data extraction
fileIDsEY =dir(strcat(pathEY,'*.csv'));
fileIDsEY = struct2cell(fileIDsEY)';
fileIDsEY(:,[2:6]) = [];

fileIDsMA =dir(strcat(pathMA,'*.csv'));
fileIDsMA = struct2cell(fileIDsMA)';
fileIDsMA(:,[2:6]) = [];

fileIDsTW = [dir(strcat(pathTW,'*.csv')); dir(strcat(pathTW,'*.xlsx'))];
fileIDsTW = struct2cell(fileIDsTW)';
fileIDsTW(:,[2:6]) = [];

%Predefine number of subjects for the analysis
m = 49
maxfinal = zeros(1,m);
minfinal = zeros(1,m);

% preallocate runorders for 4 category labels
RunOrderpos =zeros(999999,3,m);
RunOrderact =zeros(999999,3,m);
RunOrderfun =zeros(999999,3,m);
RunOrdertrans =zeros(999999,3,m);


%% Automatically generate and save annotation and confusion matrix in a loop

%initialise loop
c=1  %data we are looking for: 1 = xsens natural; 2 = xsens control 3 = apple
for n = 1:m   %patient number, manually insert number of patient selected here.

disp(patientnumber_XsensN(n))
Patient = char(patientnumber_XsensN(n));

idx1 = contains(fileIDsEY, patientnumber_XsensN(n));
rownum1 = find(idx1,1,'last');

idx2 = contains(fileIDsMA, patientnumber_XsensN(n));
rownum2 = find(idx2,1,'last');

idx3 = contains(fileIDsTW, patientnumber_XsensN(n));
rownum3 = find(idx3,1,'last');

%direct pathname and file directory and display raters' initials
if isempty(rownum1) == 1,
    raters = 'TW&MA'
    pn1 = pathTW;
     else pn1 = pathEY;
    raters = 'EY&MA'
end
     pn2 = pathMA;
    
%define filename in preparation for loading in following steps
if raters == 'TW&MA'
    fn1 = char(fileIDsTW(rownum3));
else fn1 = char(fileIDsEY(rownum1));
end
     fn2 = char(fileIDsMA(rownum2));
     
%break here!!!!!!!!!!!!!!!!!!
%load rater 1's rating
F1 = pwd;

cd(pn1)
fprintf('CHOSEN FOLDER No.1: %s\n', pn1)
fprintf('\tReading: ');
rawClassData1 = readtable(fn1);
fprintf('%s [%d x %d] & ', fn1);

%load rater 2's rating
cd(pn2)
fprintf('\nCHOSEN FOLDER No.2: %s\n', pn2)
fprintf('\tReading: ');
rawClassData2 = readtable(fn2);
fprintf('%s [%d x %d] & ', fn2);

cd(F1)

clear pathEY pathMA pathTW idx1 idx2 idx3 pn1 pn2 rownum1 rownum2 rownum3 fileIDs*

% data preparation    BREAK HERE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%rater1
ClassData1 = rawClassData1(:,[7 11 12]);
ClassData1 =table2array(ClassData1);
%rater2
ClassData2 = rawClassData2(:,[7 11 12]);
ClassData2 =table2array(ClassData2);

%Convert Rater 1 annotation seconds into DeciSeconds (this is to preserve those classifications that start and end within same second)
ClassDataDS1=[ClassData1(:,[2 3])*10 ClassData1(:,1) zeros(length(ClassData1),1)];

%Convert Rater 2 annotation seconds into DeciSeconds (this is to preserve those classifications that start and end within same second)
ClassDataDS2=[ClassData2(:,[2 3])*10 ClassData2(:,1) zeros(length(ClassData2),1)];

clear rawClassData1 rawClassData2 ClassData1 ClassData2 pn1 pn2 F1

%Correct for NPT unsync
%ClassDataDS2(:,1:2) = ClassDataDS2(:,1:2)-58;

% Generate Running Order of Annotations 
ClassDataDS1 =round(ClassDataDS1); %rater 1
ClassDataDS2 =round(ClassDataDS2); %rater 2

%%remove zero start index from all class files
if ClassDataDS1(1,1) ==0
   ClassDataDS1(1,1)=10;
end 
if ClassDataDS1(2,1) ==0
   ClassDataDS1(2,1)=10;
end 
if ClassDataDS2(1,1)==0
    ClassDataDS2(1,1)=10;
end 
if ClassDataDS2(2,1)==0
    ClassDataDS2(2,1)=10;
end 

% cut classification files same size
max1 =max(ClassDataDS1(:,2)); %find last annotation rater 1
max2 =max(ClassDataDS2(:,2)); %find last annotation rater 2

minfinal(1,n) = min(max1,max2); %find where to cut annotation files
maxfinal(1,n) = max(max1,max2);


% Writing into RunOrder format   %BREAKHERE!!!!!!!!!!!!!!!!
%add postures to RunOrderpos
%rater 1
for i=1:length(ClassDataDS1)
    if ClassDataDS1(i,3)<6
    RunOrderpos(ClassDataDS1(i,1): ClassDataDS1(i,2),1,n)=  ClassDataDS1(i,3);
    end
end 

%rater 2
for i=1:length(ClassDataDS2)
    if ClassDataDS2(i,3)<6
    RunOrderpos(ClassDataDS2(i,1): ClassDataDS2(i,2),2,n)=  ClassDataDS2(i,3);
    end
end 

%add functional to RunOrderfun
%rater 1
funindex1 =find(ClassDataDS1(:,3)==6 | ClassDataDS1(:,3)==8 ...
        | ClassDataDS1(:,3)==9 | ClassDataDS1(:,3)==27 | ClassDataDS1(:,3)==22 ...
        | ClassDataDS1(:,3)==23 | ClassDataDS1(:,3)==26 | ClassDataDS1(:,3)==24);
funindex1 =[1;funindex1];

for k =1:length(funindex1)-1
    t=k+1;
    RunOrderfun(ClassDataDS1(funindex1(k)):ClassDataDS1(funindex1(t)),1,n)= ClassDataDS1(funindex1(t),3);
end

%rater 2
funindex2 =find(ClassDataDS2(:,3)==6 | ClassDataDS2(:,3)==8 ...
        | ClassDataDS2(:,3)==9 | ClassDataDS2(:,3)==27 | ClassDataDS2(:,3)==22 ...
        | ClassDataDS2(:,3)==23 | ClassDataDS2(:,3)==26 | ClassDataDS2(:,3)==24);
funindex2 =[1;funindex2];
    
for k =1:length(funindex2)-1
    t=k+1;
    RunOrderfun(ClassDataDS2(funindex2(k)):ClassDataDS2(funindex2(t)),2,n)= ClassDataDS2(funindex2(t),3);
end

%
% add activities to RunOrderact
for i=1:length(ClassDataDS1)
    if ClassDataDS1(i,3)>9 && ClassDataDS1(i,3)<18
    RunOrderact(ClassDataDS1(i,1): ClassDataDS1(i,2),1,n)=  ClassDataDS1(i,3);
    end
end 

%MS
for i=1:length(ClassDataDS2)
    if ClassDataDS2(i,3)>9 && ClassDataDS2(i,3)<18
    RunOrderact(ClassDataDS2(i,1): ClassDataDS2(i,2),2,n)=  ClassDataDS2(i,3);
    end
end


%add transition to RunOrdertrans
%Rater1
for i=1:length(ClassDataDS1)
    if ClassDataDS1(i,3) ==28
    RunOrdertrans(ClassDataDS1(i,1): ClassDataDS1(i,2),1,n)=  ClassDataDS1(i,3);
    end
end

%Rater2
for i=1:length(ClassDataDS2)
    if ClassDataDS2(i,3)==28
    RunOrdertrans(ClassDataDS2(i,1): ClassDataDS2(i,2),2,n)=  ClassDataDS2(i,3);
    end
end

%BREAK HERE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% label occluded or sensor removed periods
for i=1:length(ClassDataDS1)
    if ClassDataDS1(i,3) >17 && ClassDataDS1(i,3) <21
    RunOrderact(ClassDataDS1(i,1): ClassDataDS1(i,2),1,n)= 18;  %original:= ClassDataDS1(i,3);
    end
end 


for i=1:length(ClassDataDS2)
    if ClassDataDS2(i,3) >17 && ClassDataDS2(i,3) <21
    RunOrderact(ClassDataDS2(i,1): ClassDataDS2(i,2),2,n)= 18;
    end
end

for i=1:length(ClassDataDS1)
    if ClassDataDS1(i,3) >17 && ClassDataDS1(i,3) <21
    RunOrderpos(ClassDataDS1(i,1): ClassDataDS1(i,2),1,n)= 18;
    end
end 

for i=1:length(ClassDataDS2)
    if ClassDataDS2(i,3) >17 && ClassDataDS2(i,3) <21
    RunOrderpos(ClassDataDS2(i,1): ClassDataDS2(i,2),2,n)= 18;
    end
end

for i=1:length(ClassDataDS1)
    if ClassDataDS1(i,3) >17 && ClassDataDS1(i,3) <21
    RunOrderfun(ClassDataDS1(i,1): ClassDataDS1(i,2),1,n)= 18;
    end
end 

for i=1:length(ClassDataDS2)
    if ClassDataDS2(i,3) >17 && ClassDataDS2(i,3) <21
    RunOrderfun(ClassDataDS2(i,1): ClassDataDS2(i,2),2,n)= 18;
    end
end



% Remove non-task / non-rated zero rows
% where no-task period is remain unlabeled by both raters, increasing
% agreement on the "non-labelled" categories and increase kappa
% %activities
% for i = length(RunOrderact):-1:1
%    if RunOrderact(i,1) == 0 && RunOrderact(i,2) == 0,
%        RunOrderact(i,:) = [];
%    end
% %posture
%    if RunOrderpos(i,1) == 0 && RunOrderpos(i,2) == 0,
%        RunOrderpos(i,:) = [];
%    end
% %functional
%    if RunOrderfun(i,1) == 0 && RunOrderfun(i,2) == 0,
%        RunOrderfun(i,:) = [];
%    end
% end 

% 
%Locate disagreement between two raters for plotting figures
% RunOrderpos = RunOrderpos(all(RunOrderpos(:,1:2),2),:);
RunOrderpos(:,3,n)=RunOrderpos(:,1,n)~=RunOrderpos(:,2,n);

% RunOrderfun = RunOrderfun(all(RunOrderfun(:,1:2),2),:);
RunOrderfun(:,3,n)=RunOrderfun(:,1,n)~=RunOrderfun(:,2,n);

% RunOrderact = RunOrderact(all(RunOrderact(:,1:2),2),:);
RunOrderact(:,3,n)=RunOrderact(:,1,n)~=RunOrderact(:,2,n);

% RunOrdertrans = RunOrdertrans(all(RunOrdertrans(:,1:2),2),:);
RunOrdertrans(:,3,n)=RunOrdertrans(:,1,n)~=RunOrdertrans(:,2,n);


%clear current workspace
clear ClassDataDS* fn* funindex* i max1 max2 raters Patient

%Re-initialise with path details
pathEY ='D:\Box Sync\Box Sync\Elton Shared\All task classifications (E.Y.) V3\';
pathMA ='D:\Box Sync\Box Sync\Elton Shared\All task classifications (M.A.) V3\';
pathTW ='D:\Box Sync\Box Sync\Elton Shared\Natural Xsens classifications (T.W.) V3 csvformat\';

%load list of files in each rater's folder for data extraction
fileIDsEY =dir(strcat(pathEY,'*.csv'));
fileIDsEY = struct2cell(fileIDsEY)';
fileIDsEY(:,[2:6]) = [];

fileIDsMA =dir(strcat(pathMA,'*.csv'));
fileIDsMA = struct2cell(fileIDsMA)';
fileIDsMA(:,[2:6]) = [];

fileIDsTW = [dir(strcat(pathTW,'*.csv')); dir(strcat(pathTW,'*.xlsx'))];
fileIDsTW = struct2cell(fileIDsTW)';
fileIDsTW(:,[2:6]) = [];


end


%% 
%Stagger all posture score into one matrix for a summated CM to show
%distribution and incidence of the categories
ALLRunOrderpos_XsensN = cat(1,RunOrderpos(1:minfinal(1,1),:,1),RunOrderpos(1:minfinal(1,2),:,2),...
                       RunOrderpos(1:minfinal(1,3),:,3),RunOrderpos(1:minfinal(1,4),:,4),...
                       RunOrderpos(1:minfinal(1,5),:,5),RunOrderpos(1:minfinal(1,6),:,6),...
                       RunOrderpos(1:minfinal(1,7),:,7),RunOrderpos(1:minfinal(1,8),:,8),...
                       RunOrderpos(1:minfinal(1,9),:,9),RunOrderpos(1:minfinal(1,10),:,10),...
                       RunOrderpos(1:minfinal(1,11),:,11),RunOrderpos(1:minfinal(1,12),:,12),...
                       RunOrderpos(1:minfinal(1,13),:,13),RunOrderpos(1:minfinal(1,14),:,14),...
                       RunOrderpos(1:minfinal(1,15),:,15),RunOrderpos(1:minfinal(1,16),:,16),...
                       RunOrderpos(1:minfinal(1,17),:,17),RunOrderpos(1:minfinal(1,18),:,18),...
                       RunOrderpos(1:minfinal(1,19),:,19),RunOrderpos(1:minfinal(1,20),:,20),...
                       RunOrderpos(1:minfinal(1,21),:,21),RunOrderpos(1:minfinal(1,22),:,22),...
                       RunOrderpos(1:minfinal(1,23),:,23),RunOrderpos(1:minfinal(1,24),:,24),...
                       RunOrderpos(1:minfinal(1,25),:,25),RunOrderpos(1:minfinal(1,26),:,26),...
                       RunOrderpos(1:minfinal(1,27),:,27),RunOrderpos(1:minfinal(1,28),:,28),...
                       RunOrderpos(1:minfinal(1,29),:,29),RunOrderpos(1:minfinal(1,30),:,30),...
                       RunOrderpos(1:minfinal(1,31),:,31),RunOrderpos(1:minfinal(1,32),:,32),...
                       RunOrderpos(1:minfinal(1,33),:,33),RunOrderpos(1:minfinal(1,34),:,34),...
                       RunOrderpos(1:minfinal(1,35),:,35),RunOrderpos(1:minfinal(1,36),:,36),...
                       RunOrderpos(1:minfinal(1,37),:,37),RunOrderpos(1:minfinal(1,38),:,38),...
                       RunOrderpos(1:minfinal(1,39),:,39),RunOrderpos(1:minfinal(1,40),:,40),...
                       RunOrderpos(1:minfinal(1,41),:,41),RunOrderpos(1:minfinal(1,42),:,42),...
                       RunOrderpos(1:minfinal(1,43),:,43),RunOrderpos(1:minfinal(1,44),:,44),...
                       RunOrderpos(1:minfinal(1,45),:,45),RunOrderpos(1:minfinal(1,46),:,46),...
                       RunOrderpos(1:minfinal(1,47),:,47),RunOrderpos(1:minfinal(1,48),:,48),...
                       RunOrderpos(1:minfinal(1,49),:,49));
% %                    
ALLRunOrderfun_XsensN = cat(1,RunOrderfun(1:minfinal(1,1),:,1),RunOrderfun(1:minfinal(1,2),:,2),...
                       RunOrderfun(1:minfinal(1,3),:,3),RunOrderfun(1:minfinal(1,4),:,4),...
                       RunOrderfun(1:minfinal(1,5),:,5),RunOrderfun(1:minfinal(1,6),:,6),...
                       RunOrderfun(1:minfinal(1,7),:,7),RunOrderfun(1:minfinal(1,8),:,8),...
                       RunOrderfun(1:minfinal(1,9),:,9),RunOrderfun(1:minfinal(1,10),:,10),...
                       RunOrderfun(1:minfinal(1,11),:,11),RunOrderfun(1:minfinal(1,12),:,12),...
                       RunOrderfun(1:minfinal(1,13),:,13),RunOrderfun(1:minfinal(1,14),:,14),...
                       RunOrderfun(1:minfinal(1,15),:,15),RunOrderfun(1:minfinal(1,16),:,16),...
                       RunOrderfun(1:minfinal(1,17),:,17),RunOrderfun(1:minfinal(1,18),:,18),...
                       RunOrderfun(1:minfinal(1,19),:,19),RunOrderfun(1:minfinal(1,20),:,20),...
                       RunOrderfun(1:minfinal(1,21),:,21),RunOrderfun(1:minfinal(1,22),:,22),...
                       RunOrderfun(1:minfinal(1,23),:,23),RunOrderfun(1:minfinal(1,24),:,24),...
                       RunOrderfun(1:minfinal(1,25),:,25),RunOrderfun(1:minfinal(1,26),:,26),...
                       RunOrderfun(1:minfinal(1,27),:,27),RunOrderfun(1:minfinal(1,28),:,28),...
                       RunOrderfun(1:minfinal(1,29),:,29),RunOrderfun(1:minfinal(1,30),:,30),...
                       RunOrderfun(1:minfinal(1,31),:,31),RunOrderfun(1:minfinal(1,32),:,32),...
                       RunOrderfun(1:minfinal(1,33),:,33),RunOrderfun(1:minfinal(1,34),:,34),...
                       RunOrderfun(1:minfinal(1,35),:,35),RunOrderfun(1:minfinal(1,36),:,36),...
                       RunOrderfun(1:minfinal(1,37),:,37),RunOrderfun(1:minfinal(1,38),:,38),...
                       RunOrderfun(1:minfinal(1,39),:,39),RunOrderfun(1:minfinal(1,40),:,40),...
                       RunOrderfun(1:minfinal(1,41),:,41),RunOrderfun(1:minfinal(1,42),:,42),...
                       RunOrderfun(1:minfinal(1,43),:,43),RunOrderfun(1:minfinal(1,44),:,44),...
                       RunOrderfun(1:minfinal(1,45),:,45),RunOrderfun(1:minfinal(1,46),:,46),...
                       RunOrderfun(1:minfinal(1,47),:,47),RunOrderfun(1:minfinal(1,48),:,48),...
                       RunOrderfun(1:minfinal(1,49),:,49));                 
% %                    
ALLRunOrderact_XsensN = cat(1,RunOrderact(1:minfinal(1,1),:,1),RunOrderact(1:minfinal(1,2),:,2),...
                       RunOrderact(1:minfinal(1,3),:,3),RunOrderact(1:minfinal(1,4),:,4),...
                       RunOrderact(1:minfinal(1,5),:,5),RunOrderact(1:minfinal(1,6),:,6),...
                       RunOrderact(1:minfinal(1,7),:,7),RunOrderact(1:minfinal(1,8),:,8),...
                       RunOrderact(1:minfinal(1,9),:,9),RunOrderact(1:minfinal(1,10),:,10),...
                       RunOrderact(1:minfinal(1,11),:,11),RunOrderact(1:minfinal(1,12),:,12),...
                       RunOrderact(1:minfinal(1,13),:,13),RunOrderact(1:minfinal(1,14),:,14),...
                       RunOrderact(1:minfinal(1,15),:,15),RunOrderact(1:minfinal(1,16),:,16),...
                       RunOrderact(1:minfinal(1,17),:,17),RunOrderact(1:minfinal(1,18),:,18),...
                       RunOrderact(1:minfinal(1,19),:,19),RunOrderact(1:minfinal(1,20),:,20),...
                       RunOrderact(1:minfinal(1,21),:,21),RunOrderact(1:minfinal(1,22),:,22),...
                       RunOrderact(1:minfinal(1,23),:,23),RunOrderact(1:minfinal(1,24),:,24),...
                       RunOrderact(1:minfinal(1,25),:,25),RunOrderact(1:minfinal(1,26),:,26),...
                       RunOrderact(1:minfinal(1,27),:,27),RunOrderact(1:minfinal(1,28),:,28),...
                       RunOrderact(1:minfinal(1,29),:,29),RunOrderact(1:minfinal(1,30),:,30),...
                       RunOrderact(1:minfinal(1,31),:,31),RunOrderact(1:minfinal(1,32),:,32),...
                       RunOrderact(1:minfinal(1,33),:,33),RunOrderact(1:minfinal(1,34),:,34),...
                       RunOrderact(1:minfinal(1,35),:,35),RunOrderact(1:minfinal(1,36),:,36),...
                       RunOrderact(1:minfinal(1,37),:,37),RunOrderact(1:minfinal(1,38),:,38),...
                       RunOrderact(1:minfinal(1,39),:,39),RunOrderact(1:minfinal(1,40),:,40),...
                       RunOrderact(1:minfinal(1,41),:,41),RunOrderact(1:minfinal(1,42),:,42),...
                       RunOrderact(1:minfinal(1,43),:,43),RunOrderact(1:minfinal(1,44),:,44),...
                       RunOrderact(1:minfinal(1,45),:,45),RunOrderact(1:minfinal(1,46),:,46),...
                       RunOrderact(1:minfinal(1,47),:,47),RunOrderact(1:minfinal(1,48),:,48),...
                       RunOrderact(1:minfinal(1,49),:,49));             
                   

% Remove occluded rows
ALLRunOrderact_XsensN(ALLRunOrderact_XsensN(:,1)==18,:) = [];
ALLRunOrderact_XsensN(ALLRunOrderact_XsensN(:,2)==18,:) = [];
ALLRunOrderpos_XsensN(ALLRunOrderpos_XsensN(:,1)==18,:) = [];
ALLRunOrderpos_XsensN(ALLRunOrderpos_XsensN(:,2)==18,:) = [];
ALLRunOrderfun_XsensN(ALLRunOrderfun_XsensN(:,1)==18,:) = [];
ALLRunOrderfun_XsensN(ALLRunOrderfun_XsensN(:,2)==18,:) = [];

clear RunOrder* maxfinal minfinal m n k t ClassData*
%% Repeat for XsensC 
%Predefine number of subjects for the analysis
m = 18
maxfinal = zeros(1,m);
minfinal = zeros(1,m);

% preallocate runorders for 4 category labels
RunOrderpos =zeros(999999,3,m);
RunOrderact =zeros(999999,3,m);
RunOrderfun =zeros(999999,3,m);
RunOrdertrans =zeros(999999,3,m);


%% Automatically generate and save annotation and confusion matrix in a loop - XSENSC GROUP

%initialise loop
for n = 1:m   %patient number, manually insert number of patient selected here.

disp(patientnumber_XsensC(n))
Patient = char(patientnumber_XsensC(n));

idx1 = contains(fileIDsEY, patientnumber_XsensC(n));
rownum1 = find(idx1,1,'last');

idx2 = contains(fileIDsMA, patientnumber_XsensC(n));
rownum2 = find(idx2,1,'last');

idx3 = contains(fileIDsTW, patientnumber_XsensC(n));
rownum3 = find(idx3,1,'last');

%direct pathname and file directory and display raters' initials
if isempty(rownum1) == 1,
    raters = 'TW&MA'
    pn1 = pathTW;
     else pn1 = pathEY;
    raters = 'EY&MA'
end
     pn2 = pathMA;
    
%define filename in preparation for loading in following steps
if raters == 'TW&MA'
    fn1 = char(fileIDsTW(rownum3));
else fn1 = char(fileIDsEY(rownum1));
end
     fn2 = char(fileIDsMA(rownum2));
     
%break here!!!!!!!!!!!!!!!!!!
%load rater 1's rating
F1 = pwd;

cd(pn1)
fprintf('CHOSEN FOLDER No.1: %s\n', pn1)
fprintf('\tReading: ');
rawClassData1 = readtable(fn1);
fprintf('%s [%d x %d] & ', fn1);

%load rater 2's rating
cd(pn2)
fprintf('\nCHOSEN FOLDER No.2: %s\n', pn2)
fprintf('\tReading: ');
rawClassData2 = readtable(fn2);
fprintf('%s [%d x %d] & ', fn2);

cd(F1)

clear pathEY pathMA pathTW idx1 idx2 idx3 pn1 pn2 rownum1 rownum2 rownum3 fileIDs*

% data preparation    BREAK HERE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%rater1
ClassData1 = rawClassData1(:,[7 11 12]);
ClassData1 =table2array(ClassData1);
%rater2
ClassData2 = rawClassData2(:,[7 11 12]);
ClassData2 =table2array(ClassData2);

%Convert Rater 1 annotation seconds into DeciSeconds (this is to preserve those classifications that start and end within same second)
ClassDataDS1=[ClassData1(:,[2 3])*10 ClassData1(:,1) zeros(length(ClassData1),1)];

%Convert Rater 2 annotation seconds into DeciSeconds (this is to preserve those classifications that start and end within same second)
ClassDataDS2=[ClassData2(:,[2 3])*10 ClassData2(:,1) zeros(length(ClassData2),1)];

clear rawClassData1 rawClassData2 ClassData1 ClassData2 pn1 pn2 F1

%Correct for NPT unsync
%ClassDataDS2(:,1:2) = ClassDataDS2(:,1:2)-58;

% Generate Running Order of Annotations 
ClassDataDS1 =round(ClassDataDS1); %rater 1
ClassDataDS2 =round(ClassDataDS2); %rater 2

%%remove zero start index from all class files
if ClassDataDS1(1,1) ==0
   ClassDataDS1(1,1)=10;
end 
if ClassDataDS1(2,1) ==0
   ClassDataDS1(2,1)=10;
end 
if ClassDataDS2(1,1)==0
    ClassDataDS2(1,1)=10;
end 
if ClassDataDS2(2,1)==0
    ClassDataDS2(2,1)=10;
end 

%%remove negative start index from all class files
if ClassDataDS1(1,1) <0
   ClassDataDS1(1,1)=10;
end 
if ClassDataDS1(2,1) <0
   ClassDataDS1(2,1)=10;
end 
if ClassDataDS2(1,1) <0
    ClassDataDS2(1,1)=10;
end 
if ClassDataDS2(2,1) <0
    ClassDataDS2(2,1)=10;
end 

% cut classification files same size
max1 =max(ClassDataDS1(:,2)); %find last annotation rater 1
max2 =max(ClassDataDS2(:,2)); %find last annotation rater 2

minfinal(1,n) = min(max1,max2); %find where to cut annotation files
maxfinal(1,n) = max(max1,max2);


% Writing into RunOrder format   %BREAKHERE!!!!!!!!!!!!!!!!
%add postures to RunOrderpos
%rater 1
for i=1:length(ClassDataDS1)
    if ClassDataDS1(i,3)<6
    RunOrderpos(ClassDataDS1(i,1): ClassDataDS1(i,2),1,n)=  ClassDataDS1(i,3);
    end
end 

%rater 2
for i=1:length(ClassDataDS2)
    if ClassDataDS2(i,3)<6
    RunOrderpos(ClassDataDS2(i,1): ClassDataDS2(i,2),2,n)=  ClassDataDS2(i,3);
    end
end 

%add functional to RunOrderfun
%rater 1
funindex1 =find(ClassDataDS1(:,3)==6 | ClassDataDS1(:,3)==8 ...
        | ClassDataDS1(:,3)==9 | ClassDataDS1(:,3)==27 | ClassDataDS1(:,3)==22 ...
        | ClassDataDS1(:,3)==23 | ClassDataDS1(:,3)==26 | ClassDataDS1(:,3)==24);
funindex1 =[1;funindex1];

for k =1:length(funindex1)-1
    t=k+1;
    RunOrderfun(ClassDataDS1(funindex1(k)):ClassDataDS1(funindex1(t)),1,n)= ClassDataDS1(funindex1(t),3);
end

%rater 2
funindex2 =find(ClassDataDS2(:,3)==6 | ClassDataDS2(:,3)==8 ...
        | ClassDataDS2(:,3)==9 | ClassDataDS2(:,3)==27 | ClassDataDS2(:,3)==22 ...
        | ClassDataDS2(:,3)==23 | ClassDataDS2(:,3)==26 | ClassDataDS2(:,3)==24);
funindex2 =[1;funindex2];
    
for k =1:length(funindex2)-1
    t=k+1;
    RunOrderfun(ClassDataDS2(funindex2(k)):ClassDataDS2(funindex2(t)),2,n)= ClassDataDS2(funindex2(t),3);
end

%
% add activities to RunOrderact
for i=1:length(ClassDataDS1)
    if ClassDataDS1(i,3)>9 && ClassDataDS1(i,3)<18
    RunOrderact(ClassDataDS1(i,1): ClassDataDS1(i,2),1,n)=  ClassDataDS1(i,3);
    end
end 

%MS
for i=1:length(ClassDataDS2)
    if ClassDataDS2(i,3)>9 && ClassDataDS2(i,3)<18
    RunOrderact(ClassDataDS2(i,1): ClassDataDS2(i,2),2,n)=  ClassDataDS2(i,3);
    end
end


%add transition to RunOrdertrans
%Rater1
for i=1:length(ClassDataDS1)
    if ClassDataDS1(i,3) ==28
    RunOrdertrans(ClassDataDS1(i,1): ClassDataDS1(i,2),1,n)=  ClassDataDS1(i,3);
    end
end

%Rater2
for i=1:length(ClassDataDS2)
    if ClassDataDS2(i,3)==28
    RunOrdertrans(ClassDataDS2(i,1): ClassDataDS2(i,2),2,n)=  ClassDataDS2(i,3);
    end
end

%BREAK HERE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% label occluded or sensor removed periods
for i=1:length(ClassDataDS1)
    if ClassDataDS1(i,3) >17 && ClassDataDS1(i,3) <21
    RunOrderact(ClassDataDS1(i,1): ClassDataDS1(i,2),1,n)= 18;  %original:= ClassDataDS1(i,3);
    end
end 


for i=1:length(ClassDataDS2)
    if ClassDataDS2(i,3) >17 && ClassDataDS2(i,3) <21
    RunOrderact(ClassDataDS2(i,1): ClassDataDS2(i,2),2,n)= 18;
    end
end

for i=1:length(ClassDataDS1)
    if ClassDataDS1(i,3) >17 && ClassDataDS1(i,3) <21
    RunOrderpos(ClassDataDS1(i,1): ClassDataDS1(i,2),1,n)= 18;
    end
end 

for i=1:length(ClassDataDS2)
    if ClassDataDS2(i,3) >17 && ClassDataDS2(i,3) <21
    RunOrderpos(ClassDataDS2(i,1): ClassDataDS2(i,2),2,n)= 18;
    end
end

for i=1:length(ClassDataDS1)
    if ClassDataDS1(i,3) >17 && ClassDataDS1(i,3) <21
    RunOrderfun(ClassDataDS1(i,1): ClassDataDS1(i,2),1,n)= 18;
    end
end 

for i=1:length(ClassDataDS2)
    if ClassDataDS2(i,3) >17 && ClassDataDS2(i,3) <21
    RunOrderfun(ClassDataDS2(i,1): ClassDataDS2(i,2),2,n)= 18;
    end
end


% 
%Locate disagreement between two raters for plotting figures
% RunOrderpos = RunOrderpos(all(RunOrderpos(:,1:2),2),:);
RunOrderpos(:,3,n)=RunOrderpos(:,1,n)~=RunOrderpos(:,2,n);

% RunOrderfun = RunOrderfun(all(RunOrderfun(:,1:2),2),:);
RunOrderfun(:,3,n)=RunOrderfun(:,1,n)~=RunOrderfun(:,2,n);

% RunOrderact = RunOrderact(all(RunOrderact(:,1:2),2),:);
RunOrderact(:,3,n)=RunOrderact(:,1,n)~=RunOrderact(:,2,n);

% RunOrdertrans = RunOrdertrans(all(RunOrdertrans(:,1:2),2),:);
RunOrdertrans(:,3,n)=RunOrdertrans(:,1,n)~=RunOrdertrans(:,2,n);


%clear current workspace
clear ClassDataDS* fn* funindex* i max1 max2 raters Patient

%Re-initialise with path details
pathEY ='D:\Box Sync\Box Sync\Elton Shared\All task classifications (E.Y.) V3\';
pathMA ='D:\Box Sync\Box Sync\Elton Shared\All task classifications (M.A.) V3\';
pathTW ='D:\Box Sync\Box Sync\Elton Shared\Natural Xsens classifications (T.W.) V3 csvformat\';

%load list of files in each rater's folder for data extraction
fileIDsEY =dir(strcat(pathEY,'*.csv'));
fileIDsEY = struct2cell(fileIDsEY)';
fileIDsEY(:,[2:6]) = [];

fileIDsMA =dir(strcat(pathMA,'*.csv'));
fileIDsMA = struct2cell(fileIDsMA)';
fileIDsMA(:,[2:6]) = [];

fileIDsTW = [dir(strcat(pathTW,'*.csv')); dir(strcat(pathTW,'*.xlsx'))];
fileIDsTW = struct2cell(fileIDsTW)';
fileIDsTW(:,[2:6]) = [];


end


%% 
%Stagger all posture score into one matrix for a summated CM to show
%distribution and incidence of the categories
ALLRunOrderpos_XsensC = cat(1,RunOrderpos(1:minfinal(1,1),:,1),RunOrderpos(1:minfinal(1,2),:,2),...
                       RunOrderpos(1:minfinal(1,3),:,3),RunOrderpos(1:minfinal(1,4),:,4),...
                       RunOrderpos(1:minfinal(1,5),:,5),RunOrderpos(1:minfinal(1,6),:,6),...
                       RunOrderpos(1:minfinal(1,7),:,7),RunOrderpos(1:minfinal(1,8),:,8),...
                       RunOrderpos(1:minfinal(1,9),:,9),RunOrderpos(1:minfinal(1,10),:,10),...
                       RunOrderpos(1:minfinal(1,11),:,11),RunOrderpos(1:minfinal(1,12),:,12),...
                       RunOrderpos(1:minfinal(1,13),:,13),RunOrderpos(1:minfinal(1,14),:,14),...
                       RunOrderpos(1:minfinal(1,15),:,15),RunOrderpos(1:minfinal(1,16),:,16),...
                       RunOrderpos(1:minfinal(1,17),:,17),RunOrderpos(1:minfinal(1,18),:,18));

% %                    
ALLRunOrderfun_XsensC = cat(1,RunOrderfun(1:minfinal(1,1),:,1),RunOrderfun(1:minfinal(1,2),:,2),...
                       RunOrderfun(1:minfinal(1,3),:,3),RunOrderfun(1:minfinal(1,4),:,4),...
                       RunOrderfun(1:minfinal(1,5),:,5),RunOrderfun(1:minfinal(1,6),:,6),...
                       RunOrderfun(1:minfinal(1,7),:,7),RunOrderfun(1:minfinal(1,8),:,8),...
                       RunOrderfun(1:minfinal(1,9),:,9),RunOrderfun(1:minfinal(1,10),:,10),...
                       RunOrderfun(1:minfinal(1,11),:,11),RunOrderfun(1:minfinal(1,12),:,12),...
                       RunOrderfun(1:minfinal(1,13),:,13),RunOrderfun(1:minfinal(1,14),:,14),...
                       RunOrderfun(1:minfinal(1,15),:,15),RunOrderfun(1:minfinal(1,16),:,16),...
                       RunOrderfun(1:minfinal(1,17),:,17),RunOrderfun(1:minfinal(1,18),:,18));
                   
% %                    
ALLRunOrderact_XsensC = cat(1,RunOrderact(1:minfinal(1,1),:,1),RunOrderact(1:minfinal(1,2),:,2),...
                       RunOrderact(1:minfinal(1,3),:,3),RunOrderact(1:minfinal(1,4),:,4),...
                       RunOrderact(1:minfinal(1,5),:,5),RunOrderact(1:minfinal(1,6),:,6),...
                       RunOrderact(1:minfinal(1,7),:,7),RunOrderact(1:minfinal(1,8),:,8),...
                       RunOrderact(1:minfinal(1,9),:,9),RunOrderact(1:minfinal(1,10),:,10),...
                       RunOrderact(1:minfinal(1,11),:,11),RunOrderact(1:minfinal(1,12),:,12),...
                       RunOrderact(1:minfinal(1,13),:,13),RunOrderact(1:minfinal(1,14),:,14),...
                       RunOrderact(1:minfinal(1,15),:,15),RunOrderact(1:minfinal(1,16),:,16),...
                       RunOrderact(1:minfinal(1,17),:,17),RunOrderact(1:minfinal(1,18),:,18));
                      

% Remove occluded rows
ALLRunOrderact_XsensC(ALLRunOrderact_XsensC(:,1)==18,:) = [];
ALLRunOrderact_XsensC(ALLRunOrderact_XsensC(:,2)==18,:) = [];
ALLRunOrderpos_XsensC(ALLRunOrderpos_XsensC(:,1)==18,:) = [];
ALLRunOrderpos_XsensC(ALLRunOrderpos_XsensC(:,2)==18,:) = [];
ALLRunOrderfun_XsensC(ALLRunOrderfun_XsensC(:,1)==18,:) = [];
ALLRunOrderfun_XsensC(ALLRunOrderfun_XsensC(:,2)==18,:) = [];


clear RunOrder* k m t n maxfinal minfinal m maxfinal minfinal 

%% Repeat for Apple 
%Predefine number of subjects for the analysis
m = 9
maxfinal = zeros(1,m);
minfinal = zeros(1,m);

% preallocate runorders for 4 category labels
RunOrderpos =zeros(999999,3,m);
RunOrderact =zeros(999999,3,m);
RunOrderfun =zeros(999999,3,m);
RunOrdertrans =zeros(999999,3,m);

%Change MA pathname for apple
pathMA ='D:\Box Sync\Box Sync\Elton Shared\All task classifications (M.A.) V3\Renamed Apple\';
%Redirect to new apple MA folder 
fileIDsMA =dir(strcat(pathMA,'*.csv'));
fileIDsMA = struct2cell(fileIDsMA)';
fileIDsMA(:,[2:6]) = [];


%% Automatically generate and save annotation and confusion matrix in a loop  - APPLE GROUP

%initialise loop
c=1  %data we are looking for: 1 = xsens natural; 2 = xsens control 3 = apple
for n = 1:m   %patient number, manually insert number of patient selected here.

disp(patientnumber_Apple(n))
Patient = char(patientnumber_Apple(n));

idx1 = contains(fileIDsEY, patientnumber_Apple(n));
rownum1 = find(idx1,1,'last');

idx2 = contains(fileIDsMA, patientnumber_Apple(n));
rownum2 = find(idx2,1,'last');

idx3 = contains(fileIDsTW, patientnumber_Apple(n));
rownum3 = find(idx3,1,'last');

%direct pathname and file directory and display raters' initials
if isempty(rownum1) == 1,
    raters = 'TW&MA'
    pn1 = pathTW;
     else pn1 = pathEY;
    raters = 'EY&MA'
end
     pn2 = pathMA;
    
%define filename in preparation for loading in following steps
if raters == 'TW&MA'
    fn1 = char(fileIDsTW(rownum3));
else fn1 = char(fileIDsEY(rownum1));
end
     fn2 = char(fileIDsMA(rownum2));
     
%break here!!!!!!!!!!!!!!!!!!
%load rater 1's rating
F1 = pwd;

cd(pn1)
fprintf('CHOSEN FOLDER No.1: %s\n', pn1)
fprintf('\tReading: ');
rawClassData1 = readtable(fn1);
fprintf('%s [%d x %d] & ', fn1);

%load rater 2's rating
cd(pn2)
fprintf('\nCHOSEN FOLDER No.2: %s\n', pn2)
fprintf('\tReading: ');
rawClassData2 = readtable(fn2);
fprintf('%s [%d x %d] & ', fn2);

cd(F1)

clear pathEY pathMA pathTW idx1 idx2 idx3 pn1 pn2 rownum1 rownum2 rownum3 fileIDs*

% data preparation    BREAK HERE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%rater1
ClassData1 = rawClassData1(:,[7 11 12]);
ClassData1 =table2array(ClassData1);
%rater2
ClassData2 = rawClassData2(:,[7 11 12]);
ClassData2 =table2array(ClassData2);

%Convert Rater 1 annotation seconds into DeciSeconds (this is to preserve those classifications that start and end within same second)
ClassDataDS1=[ClassData1(:,[2 3])*10 ClassData1(:,1) zeros(length(ClassData1),1)];

%Convert Rater 2 annotation seconds into DeciSeconds (this is to preserve those classifications that start and end within same second)
ClassDataDS2=[ClassData2(:,[2 3])*10 ClassData2(:,1) zeros(length(ClassData2),1)];

clear rawClassData1 rawClassData2 ClassData1 ClassData2 pn1 pn2 F1

%Correct for NPT unsync
%ClassDataDS2(:,1:2) = ClassDataDS2(:,1:2)-58;

% Generate Running Order of Annotations 
ClassDataDS1 =round(ClassDataDS1); %rater 1
ClassDataDS2 =round(ClassDataDS2); %rater 2

%%remove zero start index from all class files
if ClassDataDS1(1,1) ==0
   ClassDataDS1(1,1)=10;
end 
if ClassDataDS1(2,1) ==0
   ClassDataDS1(2,1)=10;
end 
if ClassDataDS2(1,1)==0
    ClassDataDS2(1,1)=10;
end 
if ClassDataDS2(2,1)==0
    ClassDataDS2(2,1)=10;
end 

% cut classification files same size
max1 =max(ClassDataDS1(:,2)); %find last annotation rater 1
max2 =max(ClassDataDS2(:,2)); %find last annotation rater 2

minfinal(1,n) = min(max1,max2); %find where to cut annotation files
maxfinal(1,n) = max(max1,max2);


% Writing into RunOrder format   %BREAKHERE!!!!!!!!!!!!!!!!
%add postures to RunOrderpos
%rater 1
for i=1:length(ClassDataDS1)
    if ClassDataDS1(i,3)<6
    RunOrderpos(ClassDataDS1(i,1): ClassDataDS1(i,2),1,n)=  ClassDataDS1(i,3);
    end
end 

%rater 2
for i=1:length(ClassDataDS2)
    if ClassDataDS2(i,3)<6
    RunOrderpos(ClassDataDS2(i,1): ClassDataDS2(i,2),2,n)=  ClassDataDS2(i,3);
    end
end 

%add functional to RunOrderfun
%rater 1
funindex1 =find(ClassDataDS1(:,3)==6 | ClassDataDS1(:,3)==8 ...
        | ClassDataDS1(:,3)==9 | ClassDataDS1(:,3)==27 | ClassDataDS1(:,3)==22 ...
        | ClassDataDS1(:,3)==23 | ClassDataDS1(:,3)==26 | ClassDataDS1(:,3)==24);
funindex1 =[1;funindex1];

for k =1:length(funindex1)-1
    t=k+1;
    RunOrderfun(ClassDataDS1(funindex1(k)):ClassDataDS1(funindex1(t)),1,n)= ClassDataDS1(funindex1(t),3);
end

%rater 2
funindex2 =find(ClassDataDS2(:,3)==6 | ClassDataDS2(:,3)==8 ...
        | ClassDataDS2(:,3)==9 | ClassDataDS2(:,3)==27 | ClassDataDS2(:,3)==22 ...
        | ClassDataDS2(:,3)==23 | ClassDataDS2(:,3)==26 | ClassDataDS2(:,3)==24);
funindex2 =[1;funindex2];
    
for k =1:length(funindex2)-1
    t=k+1;
    RunOrderfun(ClassDataDS2(funindex2(k)):ClassDataDS2(funindex2(t)),2,n)= ClassDataDS2(funindex2(t),3);
end

%
% add activities to RunOrderact
for i=1:length(ClassDataDS1)
    if ClassDataDS1(i,3)>9 && ClassDataDS1(i,3)<18
    RunOrderact(ClassDataDS1(i,1): ClassDataDS1(i,2),1,n)=  ClassDataDS1(i,3);
    end
end 

%MS
for i=1:length(ClassDataDS2)
    if ClassDataDS2(i,3)>9 && ClassDataDS2(i,3)<18
    RunOrderact(ClassDataDS2(i,1): ClassDataDS2(i,2),2,n)=  ClassDataDS2(i,3);
    end
end


%add transition to RunOrdertrans
%Rater1
for i=1:length(ClassDataDS1)
    if ClassDataDS1(i,3) ==28
    RunOrdertrans(ClassDataDS1(i,1): ClassDataDS1(i,2),1,n)=  ClassDataDS1(i,3);
    end
end

%Rater2
for i=1:length(ClassDataDS2)
    if ClassDataDS2(i,3)==28
    RunOrdertrans(ClassDataDS2(i,1): ClassDataDS2(i,2),2,n)=  ClassDataDS2(i,3);
    end
end

%BREAK HERE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% label occluded or sensor removed periods
for i=1:length(ClassDataDS1)
    if ClassDataDS1(i,3) >17 && ClassDataDS1(i,3) <21
    RunOrderact(ClassDataDS1(i,1): ClassDataDS1(i,2),1,n)= 18;  %original:= ClassDataDS1(i,3);
    end
end 


for i=1:length(ClassDataDS2)
    if ClassDataDS2(i,3) >17 && ClassDataDS2(i,3) <21
    RunOrderact(ClassDataDS2(i,1): ClassDataDS2(i,2),2,n)= 18;
    end
end

for i=1:length(ClassDataDS1)
    if ClassDataDS1(i,3) >17 && ClassDataDS1(i,3) <21
    RunOrderpos(ClassDataDS1(i,1): ClassDataDS1(i,2),1,n)= 18;
    end
end 

for i=1:length(ClassDataDS2)
    if ClassDataDS2(i,3) >17 && ClassDataDS2(i,3) <21
    RunOrderpos(ClassDataDS2(i,1): ClassDataDS2(i,2),2,n)= 18;
    end
end

for i=1:length(ClassDataDS1)
    if ClassDataDS1(i,3) >17 && ClassDataDS1(i,3) <21
    RunOrderfun(ClassDataDS1(i,1): ClassDataDS1(i,2),1,n)= 18;
    end
end 

for i=1:length(ClassDataDS2)
    if ClassDataDS2(i,3) >17 && ClassDataDS2(i,3) <21
    RunOrderfun(ClassDataDS2(i,1): ClassDataDS2(i,2),2,n)= 18;
    end
end


% 
%Locate disagreement between two raters for plotting figures
% RunOrderpos = RunOrderpos(all(RunOrderpos(:,1:2),2),:);
RunOrderpos(:,3,n)=RunOrderpos(:,1,n)~=RunOrderpos(:,2,n);

% RunOrderfun = RunOrderfun(all(RunOrderfun(:,1:2),2),:);
RunOrderfun(:,3,n)=RunOrderfun(:,1,n)~=RunOrderfun(:,2,n);

% RunOrderact = RunOrderact(all(RunOrderact(:,1:2),2),:);
RunOrderact(:,3,n)=RunOrderact(:,1,n)~=RunOrderact(:,2,n);

% RunOrdertrans = RunOrdertrans(all(RunOrdertrans(:,1:2),2),:);
RunOrdertrans(:,3,n)=RunOrdertrans(:,1,n)~=RunOrdertrans(:,2,n);


%clear current workspace
clear ClassDataDS* fn* funindex* i max1 max2 raters Patient

%Re-initialise with path details
pathEY ='D:\Box Sync\Box Sync\Elton Shared\All task classifications (E.Y.) V3\';
pathMA ='D:\Box Sync\Box Sync\Elton Shared\All task classifications (M.A.) V3\Renamed Apple\';
pathTW ='D:\Box Sync\Box Sync\Elton Shared\Natural Xsens classifications (T.W.) V3 csvformat\';

%load list of files in each rater's folder for data extraction
fileIDsEY =dir(strcat(pathEY,'*.csv'));
fileIDsEY = struct2cell(fileIDsEY)';
fileIDsEY(:,[2:6]) = [];

fileIDsMA =dir(strcat(pathMA,'*.csv'));
fileIDsMA = struct2cell(fileIDsMA)';
fileIDsMA(:,[2:6]) = [];

fileIDsTW = [dir(strcat(pathTW,'*.csv')); dir(strcat(pathTW,'*.xlsx'))];
fileIDsTW = struct2cell(fileIDsTW)';
fileIDsTW(:,[2:6]) = [];


end

clear RunOrder* m k n maxfinal minfinal c patientnumber* t
%% 
%Stagger all posture score into one matrix for a summated CM to show
%distribution and incidence of the categories
ALLRunOrderpos_Apple = cat(1,RunOrderpos(1:minfinal(1,1),:,1),RunOrderpos(1:minfinal(1,2),:,2),...
                       RunOrderpos(1:minfinal(1,3),:,3),RunOrderpos(1:minfinal(1,4),:,4),...
                       RunOrderpos(1:minfinal(1,5),:,5),RunOrderpos(1:minfinal(1,6),:,6),...
                       RunOrderpos(1:minfinal(1,7),:,7),RunOrderpos(1:minfinal(1,8),:,8),...
                       RunOrderpos(1:minfinal(1,9),:,9));

% %                    
ALLRunOrderfun_Apple = cat(1,RunOrderfun(1:minfinal(1,1),:,1),RunOrderfun(1:minfinal(1,2),:,2),...
                       RunOrderfun(1:minfinal(1,3),:,3),RunOrderfun(1:minfinal(1,4),:,4),...
                       RunOrderfun(1:minfinal(1,5),:,5),RunOrderfun(1:minfinal(1,6),:,6),...
                       RunOrderfun(1:minfinal(1,7),:,7),RunOrderfun(1:minfinal(1,8),:,8),...
                       RunOrderfun(1:minfinal(1,9),:,9));
% %                    
ALLRunOrderact_Apple = cat(1,RunOrderact(1:minfinal(1,1),:,1),RunOrderact(1:minfinal(1,2),:,2),...
                       RunOrderact(1:minfinal(1,3),:,3),RunOrderact(1:minfinal(1,4),:,4),...
                       RunOrderact(1:minfinal(1,5),:,5),RunOrderact(1:minfinal(1,6),:,6),...
                       RunOrderact(1:minfinal(1,7),:,7),RunOrderact(1:minfinal(1,8),:,8),...
                       RunOrderact(1:minfinal(1,9),:,9));
                      

% Remove occluded rows
ALLRunOrderact_Apple(ALLRunOrderact_Apple(:,1)==18,:) = [];
ALLRunOrderact_Apple(ALLRunOrderact_Apple(:,2)==18,:) = [];
ALLRunOrderpos_Apple(ALLRunOrderpos_Apple(:,1)==18,:) = [];
ALLRunOrderpos_Apple(ALLRunOrderpos_Apple(:,2)==18,:) = [];
ALLRunOrderfun_Apple(ALLRunOrderfun_Apple(:,1)==18,:) = [];
ALLRunOrderfun_Apple(ALLRunOrderfun_Apple(:,2)==18,:) = [];



%%

%generate confusion matrix for posture, functional, and activities
figure(1)
ALLRunOrderpos_XsensN(ALLRunOrderpos_XsensN ==0)=99;  %label non-rated packages as 99 for cm
cm =confusionchart(ALLRunOrderpos_XsensN(:,1),ALLRunOrderpos_XsensN(:,2));
cm.Normalization = 'total-normalized';
cm.XLabel = 'Rater 1';
cm.YLabel = 'Rater 2';
cm.Title =  'Xsens Natural Summated CM Postures';
cm.RowSummary = 'row-normalized';
cm.ColumnSummary = 'column-normalized';
cmpos_XsensN = cm.NormalizedValues;
%
figure(2)
ALLRunOrderfun_XsensN(ALLRunOrderfun_XsensN ==0)=99;  %label non-rated packages as 99 for cm
cm =confusionchart(ALLRunOrderfun_XsensN(:,1),ALLRunOrderfun_XsensN(:,2));
cm.Normalization = 'total-normalized';
cm.XLabel = 'Rater 1';
cm.YLabel = 'Rater 2';
cm.Title = 'Xsens Natural Summated CM Functional';
cm.RowSummary = 'row-normalized';
cm.ColumnSummary = 'column-normalized';
cmfun_XsensN = cm.NormalizedValues;

figure(3)
ALLRunOrderact_XsensN(ALLRunOrderact_XsensN ==0)=99;  %label non-rated packages as 99 for cm
cm =confusionchart(ALLRunOrderact_XsensN(:,1),ALLRunOrderact_XsensN(:,2));
cm.Normalization = 'total-normalized';
cm.XLabel = 'Rater 1';
cm.YLabel = 'Rater 2';
cm.Title = 'Xsens Natural Summated CM Activities';
cm.RowSummary = 'row-normalized';
cm.ColumnSummary = 'column-normalized';
cmact_XsensN = cm.NormalizedValues;
% 

figure(4)
ALLRunOrderpos_XsensC(ALLRunOrderpos_XsensC ==0)=99;  %label non-rated packages as 99 for cm
cm =confusionchart(ALLRunOrderpos_XsensC(:,1),ALLRunOrderpos_XsensC(:,2));
cm.Normalization = 'total-normalized';
cm.XLabel = 'Rater 1';
cm.YLabel = 'Rater 2';
cm.Title =  'Xsens Controlled Summated CM Postures';
cm.RowSummary = 'row-normalized';
cm.ColumnSummary = 'column-normalized';
cmpos_XsensC = cm.NormalizedValues;

figure(5)
ALLRunOrderfun_XsensC(ALLRunOrderfun_XsensC ==0)=99;  %label non-rated packages as 99 for cm
cm =confusionchart(ALLRunOrderfun_XsensC(:,1),ALLRunOrderfun_XsensC(:,2));
cm.Normalization = 'total-normalized';
cm.XLabel = 'Rater 1';
cm.YLabel = 'Rater 2';
cm.Title = 'Xsens Controlled Summated CM Functional';
cm.RowSummary = 'row-normalized';
cm.ColumnSummary = 'column-normalized';
cmfun_XsensC = cm.NormalizedValues;

figure(6)
ALLRunOrderact_XsensC(ALLRunOrderact_XsensC ==0)=99;  %label non-rated packages as 99 for cm
cm =confusionchart(ALLRunOrderact_XsensC(:,1),ALLRunOrderact_XsensC(:,2));
cm.Normalization = 'total-normalized';
cm.XLabel = 'Rater 1';
cm.YLabel = 'Rater 2';
cm.Title = 'Xsens Controlled Summated CM Activities';
cm.RowSummary = 'row-normalized';
cm.ColumnSummary = 'column-normalized';
cmact_XsensC = cm.NormalizedValues;

figure(7)
ALLRunOrderpos_Apple(ALLRunOrderpos_Apple ==0)=99;  %label non-rated packages as 99 for cm
cm =confusionchart(ALLRunOrderpos_Apple(:,1),ALLRunOrderpos_Apple(:,2));
cm.Normalization = 'total-normalized';
cm.XLabel = 'Rater 1';
cm.YLabel = 'Rater 2';
cm.Title =  'Apple Summated CM Postures';
cm.RowSummary = 'row-normalized';
cm.ColumnSummary = 'column-normalized';
cmpos_Apple = cm.NormalizedValues;

figure(8)
ALLRunOrderfun_Apple(ALLRunOrderfun_Apple==0)=99;  %label non-rated packages as 99 for cm
cm =confusionchart(ALLRunOrderfun_Apple(:,1),ALLRunOrderfun_Apple(:,2));
cm.Normalization = 'total-normalized';
cm.XLabel = 'Rater 1';
cm.YLabel = 'Rater 2';
cm.Title = 'Apple Summated CM Functional';
cm.RowSummary = 'row-normalized';
cm.ColumnSummary = 'column-normalized';
cmfun_Apple = cm.NormalizedValues;

figure(9)
ALLRunOrderact_Apple(ALLRunOrderact_Apple ==0)=99;  %label non-rated packages as 99 for cm
cm =confusionchart(ALLRunOrderact_Apple(:,1),ALLRunOrderact_Apple(:,2));
cm.Normalization = 'total-normalized';
cm.XLabel = 'Rater 1';
cm.YLabel = 'Rater 2';
cm.Title = 'Apple Summated CM Activities';
cm.RowSummary = 'row-normalized';
cm.ColumnSummary = 'column-normalized';
cmact_Apple = cm.NormalizedValues;

clear pathEY pathMA pathTW fileIDs* cm

%% Analyse

%extract diagonal agreement score

%posture
for i = 1:size(cmpos_XsensN)
    mat_pos_score(i,1) = cmpos_XsensN(i,i);
end 

for i = 1:size(cmpos_XsensC)
    mat_pos_score(i,2) = cmpos_XsensC(i,i);
end 

for i = 1:size(cmpos_Apple)
    mat_pos_score(i,3) = cmpos_Apple(i,i);
end 

%Activity
for i = 1:size(cmact_XsensN)
    mat_act_score(i,1) = cmact_XsensN(i,i);
end 

for i = 1:size(cmact_XsensC)
    mat_act_score(i,2) = cmact_XsensC(i,i);
end 

for i = 1:size(cmact_Apple)
    mat_act_score(i,3) = cmact_Apple(i,i);
end 


%functional
for i = 1:size(cmfun_XsensN)
    mat_fun_score(i,1) = cmfun_XsensN(i,i);
end 

for i = 1:size(cmfun_XsensC)
    mat_fun_score(i,2) = cmfun_XsensC(i,i);
end 

for i = 1:size(cmfun_Apple)
    mat_fun_score(i,3) = cmfun_Apple(i,i);
end 

%clear mat_act mat_fun mat_pos 

%% Statistical comparison 
%column 1 = XsensN, 2 = XsensC, 3 = Apple

%posture 
figure(10)
subplot(3,1,1)
bar(mat_pos_score(:,1))
title('Xsens Natural')
ylim([0 0.5])
xticklabels({'1' '2' '3' '4' '5' '99'})

subplot(3,1,2)
bar(mat_pos_score(:,2),'FaceColor','r')
title('Xsens Controlled')
ylim([0 0.5])
xticklabels({'1' '2' '3' '4' '5' '99'})

subplot(3,1,3)
bar(mat_pos_score(:,3),'FaceColor','g')
title('Apple')
ylim([0 0.5])
xticklabels({'1' '2' '3' '4' '5' '99'})

%activity
figure(11)
subplot(3,1,1)
bar(mat_act_score(:,1))
title('Xsens Natural')
ylim([0 0.5])
xticklabels({'10' '11' '12' '13' '14' '15' '16' '17' '99'})

subplot(3,1,2)
bar(mat_act_score(:,2),'FaceColor','r')
title('Xsens Controlled')
ylim([0 0.5])
xticklabels({'10' '11' '12' '13' '14' '15' '16' '17' '99'})

subplot(3,1,3)
bar(mat_act_score(:,3),'FaceColor','g')
title('Apple')
ylim([0 0.5])
xticklabels({'10' '11' '12' '13' '14' '15' '16' '17' '99'})

%functional
figure(12)
subplot(3,1,1)
bar(mat_fun_score(:,1))
title('Xsens Natural')
ylim([0 0.5])
xticklabels({'9' '22' '23' '24' '27' '99'})

subplot(3,1,2)
bar(mat_fun_score(:,2),'FaceColor','r')
title('Xsens Controlled')
ylim([0 0.5])
xticklabels({'9' '22' '23' '24' '27' '99'})

subplot(3,1,3)
bar(mat_fun_score(:,3),'FaceColor','g')
title('Apple')
ylim([0 0.5])
xticklabels({'9' '22' '23' '24' '27' '99'})

%Stacked posture
figure(13)
bar(mat_pos_score(:,3),'FaceColor','g','FaceAlpha',0.5)
hold on
bar(mat_pos_score(:,1),'FaceAlpha',0.5)
hold on
bar(mat_pos_score(:,2),'FaceColor','r','FaceAlpha',0.5)
hold on
ylim([0 0.5])
xticklabels({'1' '2' '3' '4' '5' '99'})

%Stacked activity
figure(14)
bar(mat_act_score(:,3),'FaceColor','g','FaceAlpha',0.5)
hold on
bar(mat_act_score(:,1),'FaceAlpha',0.5)
hold on
bar(mat_act_score(:,2),'FaceColor','r','FaceAlpha',0.5)
hold on
ylim([0 0.5])
xticklabels({'10' '11' '12' '13' '14' '15' '16' '17' '99'})

%Stacked functional
figure(15)
bar(mat_fun_score(:,3),'FaceColor','g','FaceAlpha',0.5)
hold on
bar(mat_fun_score(:,1),'FaceAlpha',0.5)
hold on
bar(mat_fun_score(:,2),'FaceColor','r','FaceAlpha',0.5)
hold on
ylim([0 0.5])
xticklabels({'9' '22' '23' '24' '27' '99'})