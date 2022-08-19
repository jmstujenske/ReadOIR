function tif_name=oir2tiff(path,opts,savepath)
%
% input variables
% path: filename with path
% opts: split channels to separate tiffs,'split' or 'nosplit' (default)
% 
% output variables
% tif_name - name of tif file
%
%J.M.Stujenske, 2022
%Written based on code by Yasuhiro R. Tanaka, 2014, 2017
%
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
if nargin<2 || isempty(opts)
    opts='nosplit';
end
[folder,filename,ext]=fileparts(path);

if nargin <3||isempty(savepath)
    savepath=folder;
end
switch opts
    case 'split'
        opt=1;
    case 'nosplit'
        opt=2;
    otherwise
        error('Invalid opts type.')
end

if isempty(strfind(ext,'.oir'))
    error('input is restricted to OIR files')
end

accu_flag=0;

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
bit_depth=stdData.Metadata.bit_depth;
pixel_size=stdData.Metadata.pixel_size;
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
% n_image_tot=stdData.Metadata.sizeT_all;


%%% for old "*.oir" files
if n_ch == 0
    loc_channelId = strfind(meta,'pmt channelId');
    n_ch = numel(loc_channelId);
end
%%%
n_tz = floor(1.08e9/sizeX/sizeY/ceil(bit_depth/8)); 
if n_ch>2
    disp('Too many channels, movie with 1 or 2 channels is readable. Return blank outputs.')
%     stdData.Image{1}=[];
    return
end
        flag=0;
        [image_out,ref,index] = image_read_from_OIR(fid,sizeX,sizeY,n_tz,ref_sizeX,ref_sizeY,line_rate,flag,n_ch,accu_flag);
        if length(filelist)>2
            fid_temp=fopen(filelist(end).name,'r');fseek(fid_temp,0,'eof');flen=ftell(fid_temp);fclose(fid_temp);
            fid_temp2=fopen(filelist(end-1).name,'r');fseek(fid_temp,0,'eof');flen2=ftell(fid_temp);fclose(fid_temp);
            n_frames=index*n_ch+(index+1)*n_ch*(length(filelist)-1)+(index+1)*2-ceil((flen2-flen)/sizeX/sizeY/ceil(bit_depth/8));
        elseif length(filelist)>1
            fid_temp=fopen(filelist(end).name,'r');fseek(fid_temp,0,'eof');flen=ftell(fid_temp);fclose(fid_temp);
            n_frames=index*n_ch+(index+1)*n_ch*(length(filelist)-1)+(index+1)*2-ceil((1.077e9-flen)/sizeX/sizeY/ceil(bit_depth/8));
        else
            n_frames=index*n_ch;
        end
imagedesc=image_desc_gen(n_ch,n_z,floor(n_frames/n_ch/n_z),pixel_size,bit_depth);
        if n_ch>1
%     if ~exist(fullfile(folder,[filename,'_chan',num2str(ch_rep),'.tif']),'file')
if opt==1
    tif_name=cell(2,1);
    for ch_rep=1:n_ch
        tif_name{ch_rep}=fullfile(savepath,[filename,'_chan',num2str(ch_rep),'.tif']);
        TiffWriter{ch_rep}=Fast_BigTiff_Write(tif_name{ch_rep},pixel_size);
    end
else
    tif_name=fullfile(savepath,[filename,'.tif']);
      TiffWriter{1}=Fast_BigTiff_Write(tif_name,pixel_size,[],imagedesc);
      TiffWriter{2}=TiffWriter{1};
end
%     else
%         error('tif already exists.');
%     end

else
%     if ~exist(fullfile(folder,[filename,'.tif']),'file')
    TiffWriter{1}=Fast_BigTiff_Write(fullfile(folder,[filename,'.tif']));
%     else
%         error('tif already exists.');
%     end
        end
       for ch_rep=1:n_ch
            imimage{ch_rep}=image_out(:,:,1:index,ch_rep);
            imimage_res{ch_rep} = [];
            if n_z ~= 1
            imimage_res{ch_rep} = imimage{ch_rep}(:,:,floor(index/n_z)*n_z+1:index);
            end
        imimage{ch_rep}(:,:,floor(index/n_z)*n_z+1:index) = [];
        imimage{ch_rep} = reshape(imimage{ch_rep},sizeY,sizeX,n_z,floor(index/n_z));
         n_t=size(imimage{ch_rep},4);
         if opt==1
         for a=1:n_t;TiffWriter{ch_rep}.WriteIMG(imimage{ch_rep}(:,:,:,a)');end
         end
       end
       if opt==2
         for a=1:n_t
                for ch_rep=1:n_ch
                    TiffWriter{ch_rep}.WriteIMG(imimage{ch_rep}(:,:,:,a)');
                end
         end
       end
       fclose(fid);
        for i = 1:length(filelist)
            
            fid = fopen(filelist(i).name);
            flag=1;
            [image_out,~,index] = image_read_from_OIR(fid,sizeX,sizeY,n_tz,ref_sizeX,ref_sizeY,line_rate,flag,n_ch,accu_flag);
            for ch_rep=1:n_ch
            imimage{ch_rep} = cat(3,imimage_res{ch_rep},image_out(:,:,1:index,ch_rep));
            if ~isempty(imimage_res{ch_rep})
                index = index+size(imimage_res{ch_rep},3);
            end
            if n_z ~= 1
                imimage_res{ch_rep} = imimage{ch_rep}(:,:,floor(index/n_z)*n_z+1:index);
            end
            imimage{ch_rep}(:,:,floor(index/n_z)*n_z+1:index) = [];
            imimage{ch_rep} = reshape(imimage{ch_rep},sizeY,sizeX,n_z,floor(index/n_z));

%                 stdData_tmp = OIRxml2stdData(1,fullfile(pwd,filelist(i).name),meta,ref,imimage{1},imimage{2});
                            n_t=size(imimage{ch_rep},4);
                            if opt==1;
                for a=1:n_t;TiffWriter{ch_rep}.WriteIMG(imimage{ch_rep}(:,:,:,a)');end
                            end
            end
            if opt==2
            for a=1:n_t
                for ch_rep=1:n_ch
                    TiffWriter{ch_rep}.WriteIMG(imimage{ch_rep}(:,:,:,a)');
                end
            end
            end
        fclose(fid);
   disp(['File ',num2str(i+1),' of ',num2str(length(filelist)+1), ' converted.']);

%                 stdData.Image{2}=cat(4,stdData.Image{2},stdData_tmp.Image{2});
        end
        for ch_rep=1:n_ch
            if ~TiffWriter{ch_rep}.Closed
                close(TiffWriter{ch_rep});
            end
        end
   
%     stdData=[];
cd(oldFolder);
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

function imagedesc=image_desc_gen(n_ch,n_z,n_image_tot,pixel_size,bit_depth)
max_val=(2^bit_depth)-1;
    imagedesc=[uint8('ImageJ=1.53c'),uint8(10),...
    uint8(['images=',num2str(n_image_tot)]),uint8(10),...
    uint8(['channels=',num2str(n_ch)]),uint8(10),... 
    uint8(['slices=',num2str(n_z)]),uint8(10),... 
    uint8(['hyperstack=true']),uint8(10),... 
    uint8(['mode=grayscale']),uint8(10),... 
    uint8(['unit=micron']),uint8(10),... 
    uint8(['spacing=',num2str(pixel_size)]),uint8(10),... 
    uint8(['loop=false']),uint8(10),... 
    uint8(['min=0.0']),uint8(10),... 
    uint8(['max=',num2str(max_val),'.0']),uint8(10)];
end