% bin_mean() - Computes the mean ERP from two bins and stores it as a new  
%             bin in a GND variable. The new ERP is simply (binA+binB)/2.
%             
% Usage:
%  >> GND=bin_mean(GND,binnums,mn_bindesc,mn_ccode,mn_ccdescr);
%
% Required Inputs:
%   GND    - A GND structure variable. To create a GND variable 
%            from Kutaslab ERP files (e.g., *.mas files) use avgs2GND.m.  
%            To do the same from EEGLAB *.set files use sets2GND.m.
%            See Mass Univariate ERP Toolbox documentation for detailed  
%            information about the format of a GND variable. 
%   binnums   - vector of bins to average
%
% Optional Inputs:
%   mn_bindesc - [string] The bin descriptor for the new bin being
%                 created. {default: '(Bin #+Bin ##)/2', where # is the value
%                 binA and ## is the value of binB}
%   mn_ccode   - [integer] The condition code of the new bin being
%                 created. Condition codes are specific to Kutaslab data
%                 and can be ignored if your lab doesn't support them.
%                 {default: condition code of binA}
%   mn_ccdescr - [string] The condition code descriptor for the new bin being
%                 created. Condition codes are specific to Kutaslab data 
%                 and can be ignored if your lab doesn't support them.
%                 {default: existing descriptor for that condition code or 
%                 'Not Specified'}
%
% Notes:
% -binA and binB should index bins starting at 1 (i.e., Bin 0 is assumed to
% contain cal pulses and is stored apart from bins in GND variables)
%
% Example:
% >> GND=bin_mean(GND,2,1,'Targets+Standards'); 
%
% Author:
% David Groppe
% Kutaslab, 3/2010

%%%%%%%%%%%%%%%% REVISION LOG %%%%%%%%%%%%%%%%%
% 11/8/2016-Documentation mistakenly said this function was creating
% difference waves (rather than averaging). Thanks to Eric Fields for
% catching this.
%
% 3/15/2010-bin_dif command is added to GND variable history
%
% 12/9/2010-Original GND variable returned in case of failure (as opposed
% to NaN)

%%%%%%%%%%%%%%%% FUTURE WORK %%%%%%%%%%%%%%%%%
% - add verblevel? 

function GND=bin_mean(GND,binnums,mn_bindesc,mn_ccode,mn_ccdescr)

mn_bindesc=mn_bindesc; 

if nargin<4,
   mn_ccode=[]; 
end
if nargin<5,
   mn_ccdescr=[]; 
end

GND_copy=GND; %copy original variable in case something fails

[n_chans, n_pts, n_bins, n_subs]=size(GND.indiv_erps);
neo_bin=n_bins+1;
GND.indiv_erps(:,:,neo_bin,:)=zeros(n_chans,n_pts,1,n_subs);
use_subs=[];
for sub=1:n_subs,
            GND.indiv_erps(:,:,neo_bin,sub)=mean(GND.indiv_erps(:,:,binnums,sub),3);
            GND.indiv_bin_ct(sub,neo_bin)=-1;
            GND.indiv_bin_raw_ct(sub,neo_bin)=-1;
            use_subs=[use_subs sub];
end

    GND.sub_ct(neo_bin)=length(use_subs);
    GND.grands(:,:,neo_bin)=mean(GND.indiv_erps(:,:,neo_bin,use_subs),4);
    GND.grands_stder(:,:,neo_bin)=std(GND.indiv_erps(:,:,neo_bin,use_subs),0,4)/sqrt(GND.sub_ct(neo_bin));
    GND.grands_t(:,:,neo_bin)=GND.grands(:,:,neo_bin)./GND.grands_stder(:,:,neo_bin);
    GND.bin_info(neo_bin).bindesc=mn_bindesc;
    GND.condesc{neo_bin}=mn_ccdescr;
    GND.bin_info(neo_bin).condcode=mn_ccode;


