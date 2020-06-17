% ---------------------------------------------------------------------------
% File: compile_lumpy.m
% Purpose: Matlab script to compile the lumpy background mex files
% Notes: 
% Author:  Nick Henscheid
% Date:    9-2016, 3-2019
% Contact: nhenscheid@math.arizona.edu
% This software is in the public domain, furnished "as is", without technical
% support, and with no warranty, express or implied, as to its usefulness for
% any purpose.
% ---------------------------------------------------------------------------
disp('Compiling Lumpy Object library');

clear functions   % need to do this on windows for some reasons

addpath(genpath(pwd));

verb = 0;  % CHANGE TO 1 for DEBUGGING (VERBOSE MODE)

% Get architecture
arch=computer;
mac=strcmp(arch,'MACI64') || strcmp(arch,'MACI') || strcmp(arch,'MAC');
linux=strcmp(arch,'GLNXA64') || strcmp(arch,'GLNX86');
pc= strcmp(arch,'PCWIN64');
GPU = sign(gpuDeviceCount); % To compile with gpu or not
if(GPU)
    CUDALIB = '/usr/local/cuda/lib64';  %% CHANGE TO CUDA LIB64 DIRECTORY IF NOT HERE
    if(usejava('desktop'))
        choice = questdlg('Nvidia GPU Found!  Compile CUDA files?','Compile with CUDA?','Yes','No','Cancel Compilation','Yes');
        switch choice
            case 'Yes'
                questdlg(['Cuda library currently set to ',CUDALIB,', is this correct?'],'CUDA Library correct?','Yes','No','Cancel Compilation','Yes'); 
                GPU = 1;
            case 'No'
                GPU = 0;
            otherwise
                msgbox('Canceling compilation');
                return;
        end
    else
        GPU = input('Nvidia GPU Found!  Compile CUDA files? (1 = yes, 0 = no) ');
        input(['Compiling with CUDA.  CUDA library currently set to ',CUDALIB,', if this is incorrect modify compile_mex.m line 25! (Press any key to continue)']);
    end
end
if mac
    mexext = '.mexmaci64'; % Assuming 64 bit.
end
if linux
    mexext = '.mexa64'; % Assuming 64 bit.
end
if pc
    %mexext = '.mexw64'; 
end

%MEXOPTS = struct('CGRISRC',CGRISRC,'mexext',mexext,'CUDALIB',CUDALIB,'GPU',GPU,'VERBOSE',VERBOSE);

%GPU     = mexopts.GPU;
%CGRISRC = mexopts.CGRISRC;
%mexext  = mexopts.mexext;
%verb    = mexopts.VERBOSE;

if(GPU)
    if(verb)
        GPUFLAGS = sprintf('-v -L"%s" -lcudart -I"%s"',CUDALIB,'./source/mexsrc');
    else
        GPUFLAGS = sprintf('-L"%s" -lcudart -I"%s"',CUDALIB,'./source/mexsrc');
    end
        
end

if(verb)
    CPUFLAGS = sprintf('-v -I"%s"','./source/mexsrc');
else
    CPUFLAGS = sprintf('-I"%s"','./source/mexsrc');
end

%!!NOTE!!: Must leave a space at the beginning of each file name!

switch GPU
    case 1
        % Compile binaries for CPU/GPU
        disp('Compiling for CPU/GPU');
        str = [GPUFLAGS,' ./source/mexsrc/lumpy_mex_gpu.cu'];
        args = regexp(str,'\s+','split');
        mexcuda(args{:})    % for R2016 and above
        movefile(strcat('./lumpy_mex_gpu',mexext),strcat('./source/mexbin/lumpy_mex_gpu',mexext));
        str = [CPUFLAGS,' ./source/mexsrc/lumpy_mex_cpu.cc'];
        args = regexp(str,'\s+','split');
        mex(args{:})
        movefile(strcat('./lumpy_mex_cpu',mexext),strcat('./source/mexbin/lumpy_mex_cpu',mexext));
    case 0
        % Compile binaries for CPU only 
        disp('Compiling for CPU only');
        str = [CPUFLAGS,' ./source/mexsrc/lumpy_mex_cpu.cc'];
        if pc
            str = [CPUFLAGS,' .\source\mexsrc\lumpy_mex_cpu.cc'];
        end
        args = regexp(str,'\s+','split');
        mex(args{:})
        if pc
            movefile(strcat('.\lumpy_mex_cpu.',mexext),strcat('.\source\mexbin\lumpy_mex_cpu.',mexext));
        else
            movefile(strcat('./lumpy_mex_cpu',mexext),strcat('./source/mexbin/lumpy_mex_cpu',mexext));
        end
          
end

clear all
start