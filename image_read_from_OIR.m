function [image,ref,index] = image_read_from_OIR(fid,sizeX,sizeY,n_tz,ref_sizeX,ref_sizeY,line_rate, flag,n_ch,accu_flag)
% flag==0 -> '.oir', flag==1 -> followed files
% accu_flag,

% Copyright:
% 2014, 2017, Yasuhiro R. Tanaka; 
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
%
%Editted for speed and to tolerate large files
%Joseph M. Stujenske, 2022
%

frewind(fid)
num_line_once=ceil(30/line_rate);
num_divide=ceil(sizeY/num_line_once);
buf_all = fread(fid,inf,'uint8=>uint8');
if accu_flag~=0
    loc_0 = buf_all==0|buf_all==32;
    loc_0_s1 = [false(1,1);buf_all(1:end-1)==0|buf_all(1:end-1)==32];
    loc_0_s2 = [false(2,1);buf_all(1:end-2)==0|buf_all(1:end-2)==32];
    loc_4_s3 = [false(3,1);buf_all(1:end-3)==4];
    loc_0_s4 = [false(4,1);buf_all(1:end-4)==0|buf_all(1:end-4)==32];
    loc_95_s9 = [false(9,1);buf_all(1:end-9)==95];
    
    if num_divide>10
        loc_49_s9 = [false(9,1);buf_all(1:end-9)==49];
        loc_95_s10 = [false(10,1);buf_all(1:end-10)==95];
    end
end
fstart_p=cell(num_divide,1);

if accu_flag~=0
    for i=1:num_divide
        if i<11
            a=i;
            fstart_p{i}=find(all([loc_0,loc_0_s1,loc_0_s2,loc_4_s3,loc_0_s4,[false(8,1);buf_all(1:end-8)==47+a],loc_95_s9],2));
        else
            a=mod(i,10);
            fstart_p{i}=find(all([loc_0,loc_4_s3,[false(8,1);buf_all(1:end-8)==47+a],loc_49_s9,loc_95_s10],2));
        end
        if flag==0
            str=char(buf_all(fstart_p{i}(1)-99:fstart_p{i}(1)))';
            if contains(str,'REF')
                fstart_p{i}(1:n_ch)=[]; 
            end
        end
    end
else
        start=find(buf_all==4);
    for i=1:num_divide
        literal=['_',num2str(i-1)];
        offset=length(literal);
        mainlog=start;
        for a=0:length(literal)-1
            ins=buf_all(mainlog-offset-4+a)==uint8(literal(a+1));
            mainlog=mainlog(ins);
        end
        fstart_p{i}=mainlog+3;
        fstart_p{i}(fstart_p{i}>length(buf_all))=[];
        if flag==0
            str=char(buf_all(fstart_p{i}(1)-99:fstart_p{i}(1)))';
            if contains(str,'REF')
                fstart_p{i}(1:n_ch)=[]; 
            end
        end
    end
end
image = zeros(sizeY,sizeX,n_tz,n_ch,'uint16');
for j=1:num_divide
        for k=1:n_ch
%             end_flag=false;
           for i = 1:floor(length(fstart_p{j})/n_ch)
               end_index=fstart_p{j}(i*n_ch+k-1-(n_ch-1))+2*sizeX*min(num_line_once,sizeY-num_line_once*(j-1));
%                if end_index>length(buf_all)
%                    end_index=length(buf_all);
%                    if mod(end_index-(fstart_p{j}(i*n_ch+k-1-(n_ch-1))+1),2)==0
%                        end_index=end_index-1;
%                    end
%                    end_flag=true;
%                end
               temp=typecast(buf_all(fstart_p{j}(i*n_ch+k-1-(n_ch-1))+1:end_index),'uint16');
%                 if end_flag
%                     temp=cat(1,temp,zeros(prod([sizeX,min(num_line_once,sizeY-num_line_once*(j-1))])-length(temp),1,'uint16'));
%                 end
                image(num_line_once*(j-1)+1:min(num_line_once*j,sizeY),:,i,k) =...
                reshape(temp,[sizeX,min(num_line_once,sizeY-num_line_once*(j-1))])';
            end
        end
end

index = floor(length(fstart_p{num_divide})/n_ch);
ref=zeros(ref_sizeY,ref_sizeX,n_ch,'uint16');
end