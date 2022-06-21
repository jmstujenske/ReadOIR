function [output_list,stdData]=oir2stdData(path,mode,accu_flag)
%
% input variables
% fullpath: absolute path of file is needed. 
% all related files (sequentially recorded) is required in the same folder.
% mode == 0 -> save in separate files (< 1GB saved in -v6)
% mode == 1 -> concatenate files, CAUTION DO NOT save automatically!
% accu_flag == 0, by 1.3-fold faster than setting as 1
% accu_flag = 0 can be lead an error to read file, then use accu_flag = 1
% default setting is mode =1, accu_flag=0
%
% output variables
% output_list is a file list of saved files when mode = 0
% stdData is a struct to contain, images and metadata
% images will be found in stdData.Image{1}, etc in x*y*z*t format
%
% Copyright:
% 2014, 2017, Yasuhiro R. Tanaka; 
% License:
% This code is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published
% by the Free Software Foundation; either version 2 of the License,
% or any later version. This work is distributed in the hope that it 
% will be useful, but without any warranty; without even the implied
% warranty of merchantability or fitness for a particular purpose. See
% version 2 and version 3 of the GNU General Public License for more
% details. You should have received a copy of the GNU General Public
% License along with this program;If not, see http://www.gnu.org/licenses/.

[folder,filename,ext]=fileparts(path);

if isempty(strfind(ext,'.oir'))
    error('input is restricted to OIR files')
end

if nargin==1
    mode=1;
    accu_flag=0;
end

% search related files
if isempty(folder)
    folder = pwd;
end
oldFolder=cd(folder);

filelist = dir(pwd);
search_ph = [filename, '_[0-9]+[0-9]+[0-9]+[0-9]+[0-9]'];

for i=length(filelist):-1:1
    if isempty(regexp(filelist(i).name,search_ph,'ONCE'))
        filelist(i)=[];
    end
end

fid = fopen(path);

% metadata
[meta,~] = meta_read_from_OIR(fid);
stdData = OIRxml2stdData(0,'',meta,0);
ref_sizeX = stdData.Metadata.sizeX; % full image size (e.g. 512 etc)
ref_sizeY = stdData.Metadata.sizeY;
line_rate = stdData.Metadata.line_rate;
sizeXY= size_extract(fid,'both');
stdData.Metadata.sizeY = sizeXY(2); % ROI size
stdData.Metadata.sizeX = sizeXY(1);
n_z = stdData.Metadata.sizeZ;
sizeX = stdData.Metadata.sizeX;
sizeY = stdData.Metadata.sizeY;
loc_channelId = strfind(meta,'channel id');
n_ch = numel(loc_channelId);
%%% for old "*.oir" files
if n_ch == 0
    loc_channelId = strfind(meta,'pmt channelId');
    n_ch = numel(loc_channelId);
end
%%%
n_tz = floor(1.08e9/sizeX/sizeY/2); 
if n_ch>2
    disp('Too many channels, movie with 1 or 2 channels is readable. Return blank outputs.')
    output_list=[];
    stdData.Image{1}=[];
    return
end

switch n_ch
    case 1
        previous_index=0;
        for i = 0:1:length(filelist)
            if i==0
                flag=0;
                fullpath = path;
            else
                flag=1;
                fid = fopen(filelist(i).name);
                fullpath = fullfile(pwd,filelist(i).name);
            end
            [image1,ref,index] = image_read_from_OIR(fid,sizeX,sizeY,n_tz,ref_sizeX,ref_sizeY,line_rate,flag,n_ch,accu_flag);
            [meta,~] = meta_read_from_OIR(fid);
            image1(:,:,index+1:n_tz) = [];
            
            if i~=0
                image1 = cat(3,image_res,image1);
                if ~isempty(image_res)
                    index = index+size(image_res,3);
                end
            end
            
            image_res = [];
            if n_z ~= 1
                image_res = image1(:,:,floor(index/n_z)*n_z+1:index);% residual image will be included into next series
            end
            image1(:,:,floor(index/n_z)*n_z+1:index) = [];
            image1 = reshape(image1,sizeY,sizeX,n_z,floor(index/n_z));
            if mode==0
                % make stdData for individual files
                stdData = OIRxml2stdData(previous_index+1,fullpath,meta,ref,image1);
%                 stdData.Metadata.sizeY = size_extract(fid,'Y');
%                 stdData.Metadata.sizeX = size_extract(fid,'X');
                    sizeXY= size_extract(fid,'both');
stdData.Metadata.sizeY = sizeXY(2);% ROI size
stdData.Metadata.sizeX = sizeXY(1);
                if i<10
                    numstr = ['0', num2str(i)];
                else
                    numstr = num2str(i);
                end
                
                save(fullfile(pwd,[filename,'_',numstr,'.mat']),'stdData','-v6');
                output_list(i+1,1).name = fullfile(pwd,[filename,'_',numstr,'.mat']);
                previous_index = previous_index+size(image1,4);
                clear image1
                clear stdData
                fclose(fid);
            else
                output_list(i+1,1).name=[];
                if i==0
                    stdData = OIRxml2stdData(previous_index+1,fullpath,meta,ref,image1);
%                     stdData.Metadata.sizeY = size_extract(fid,'Y');
%                     stdData.Metadata.sizeX = size_extract(fid,'X');
                    sizeXY= size_extract(fid,'both');
                    stdData.Metadata.sizeY = sizeXY(2); % ROI size
                    stdData.Metadata.sizeX = sizeXY(1);
                else
                    stdData_tmp = OIRxml2stdData(1,fullpath,meta,ref,image1);
                    stdData.Image{1}=cat(4,stdData.Image{1},stdData_tmp.Image{1});
                end
            end
        end
        
    case 2
        output_list(1,1).name=[];
        flag=0;
        [image_out,ref,index] = image_read_from_OIR(fid,sizeX,sizeY,n_tz,ref_sizeX,ref_sizeY,line_rate,flag,n_ch,accu_flag);
       
        image1=image_out(:,:,1:index,1);
        image2=image_out(:,:,1:index,2);
        
%         image1(:,:,index+1:n_tz) = [];
%         image2(:,:,index+1:n_tz) = [];
        
        image1_res = [];
        if n_z ~= 1
            image1_res = image1(:,:,floor(index/n_z)*n_z+1:index);
        end
        
        image2_res = [];
        if n_z ~= 1
            image2_res = image2(:,:,floor(index/n_z)*n_z+1:index);
        end
        
        image1(:,:,floor(index/n_z)*n_z+1:index) = [];
        image1 = reshape(image1,sizeY,sizeX,n_z,floor(index/n_z));
        image2(:,:,floor(index/n_z)*n_z+1:index) = [];
        image2 = reshape(image2,sizeY,sizeX,n_z,floor(index/n_z));

        stdData = OIRxml2stdData(1,path,meta,ref,image1,image2);
%         stdData.Metadata.sizeY = size_extract(fid,'Y');
%         stdData.Metadata.sizeX = size_extract(fid,'X');
            sizeXY= size_extract(fid,'both');
            stdData.Metadata.sizeY = sizeXY(2); % ROI size
            stdData.Metadata.sizeX = sizeXY(1);
        n_image=cell(1,length(stdData.Image));
        for ch_rep=1:length(stdData.Image)
            n_image{ch_rep}=size(stdData.Image{ch_rep},4);
            temp=stdData.Image{ch_rep};
            stdData.Image{ch_rep}=zeros([size(stdData.Image{ch_rep},1:3),stdData.Metadata.sizeT_all],'uint16');
            stdData.Image{ch_rep}(:,:,:,1:n_image{ch_rep})=temp;
        end
        if mode ==0
            save(fullfile(pwd,[filename,'_00.mat']),'stdData','-v6');
            output_list(1,1).name = fullfile(pwd,[filename,'_00.mat']);
            previous_index = size(image1,4);
            clear image1
            clear image2
            clear stdData
            fclose(fid);
        end
     
        for i = 1:length(filelist)
            
            fid = fopen(filelist(i).name);
            flag=1;
            [image_out,~,index] = image_read_from_OIR(fid,sizeX,sizeY,n_tz,ref_sizeX,ref_sizeY,line_rate,flag,n_ch,accu_flag);
            if mode==0
                [meta,~] = meta_read_from_OIR(fid);
            end
%             image1=image_out(:,:,:,1);
%             image2=image_out(:,:,:,2);
%             image1(:,:,index+1:n_tz) = [];
            image1 = cat(3,image1_res,image_out(:,:,1:index,1));
%             image2(:,:,index+1:n_tz) = [];
            image2 = cat(3,image2_res,image_out(:,:,1:index,2));
            if ~isempty(image1_res)
                index = index+size(image1_res,3);
            end
            if n_z ~= 1
                image1_res = image1(:,:,floor(index/n_z)*n_z+1:index);
            end
            if n_z ~= 1
                image2_res = image2(:,:,floor(index/n_z)*n_z+1:index);
            end
            image1(:,:,floor(index/n_z)*n_z+1:index) = [];
            image1 = reshape(image1,sizeY,sizeX,n_z,floor(index/n_z));
            image2(:,:,floor(index/n_z)*n_z+1:index) = [];
            image2 = reshape(image2,sizeY,sizeX,n_z,floor(index/n_z));

            if mode ==0
                stdData = OIRxml2stdData(previous_index+1,fullfile(pwd,filelist(i).name),meta,ref,image1,image2);
%                 stdData.Metadata.sizeY = size_extract(fid,'Y');
%                 stdData.Metadata.sizeX = size_extract(fid,'X');
sizeXY= size_extract(fid,'both');
stdData.Metadata.sizeY = sizeXY(2); % ROI size
stdData.Metadata.sizeX = sizeXY(1);
                
                if i<10
                    numstr = ['0' num2str(i)];
                else
                    numstr = num2str(i);
                end
                
                save(fullfile(pwd,[filename,'_',numstr,'.mat']),'stdData','-v6');% no compression, fast
                output_list(i+1,1).name = fullfile(pwd,[filename,'_',numstr,'.mat']);
                previous_index = previous_index+size(image1,4);
                clear image1
                clear image2
                clear stdData
                fclose(fid);
            else
%                 stdData_tmp = OIRxml2stdData(1,fullfile(pwd,filelist(i).name),meta,ref,image1,image2);
                for ch_rep=1:2
                    switch ch_rep
                        case 1
                            n_t=size(image1,4);
                            stdData.Image{ch_rep}(:,:,:,n_image{ch_rep}+1:n_image{ch_rep}+n_t)=image1;
                        case 2
                            n_t=size(image2,4);
                            stdData.Image{ch_rep}(:,:,:,n_image{ch_rep}+1:n_image{ch_rep}+n_t)=image2;
                    end
                    n_image{ch_rep}=n_image{ch_rep}+n_t;
                end
%                 stdData.Image{2}=cat(4,stdData.Image{2},stdData_tmp.Image{2});
                output_list(i+1,1).name=[];
            end
   disp(['File ',num2str(i+1),' of ',num2str(length(filelist)+1), ' read.']);
        end
end
if mode==0
    stdData=[];
end
cd(oldFolder)
end

function sizeXY = size_extract(fid,which)
switch which
    case 'X'
        var_string{1}='width';
    case 'Y'
        var_string{1}='height';
    case 'both'
        var_string={'width','height'};
end
% sizeY
n=length(var_string);
fseek(fid,0,'bof');
loc=cell(1,n);
for a=1:5000
char_test = fread(fid,1000,'uint8=>char')';
for b=1:n
    if isempty(loc{b})
loc{b} = strfind(char_test,['<base:',var_string{b}]);
if ~isempty(loc{b})
    loc{b}=loc{b}+(1000-16)*(a-1);
end
    end
end
if ~any(cellfun('isempty',loc))
    break;
end
fseek(fid,-16,'cof');
end
% frewind(fid);
sizeXY=zeros(1,n);
for b=1:n
fseek(fid,loc{b}(1)-100,'bof');
meta_single = fread(fid,500,'uint8=>char')';
sizeXY(b) = extract_xmldata(meta_single,['base:',var_string{b}],0);
end
% sizeY
end