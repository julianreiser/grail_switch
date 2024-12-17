function [EEGnew EEGold reacmat_a reacmat_v accumat_a accumat_v trimat_a trimat_v indexstruct] = gs_reaccalc(EEG, subcondlist)

% This function calculates response times and writes everything back into
% the EEG event struct. This function also returns a matrix with
% condition-specific reaction times and accuracies.
%
% description for reacmat & accumat cols (rows always standt,walk,pert):
%  1 - single task congruent
%  2 - single task incongruent
%  3 - dual task congruent
%  4 - dual task incongruent
%  5 - dual task congruent repeat
%  6 - dual task incongruent repeat
%  7 - dual task congruent switch
%  8 - dual task incongruent switch
%
% ----- 9-16 only for reacmat -----
%
%  9 - single task congruent incorrect
% 10 - single task incongruent incorrect
% 11 - dual task congruent incorrect
% 12 - dual task incongruent incorrect
% 13 - dual task congruent repeat incorrect
% 14 - dual task incongruent repeat incorrect
% 15 - dual task congruent switch incorrect
% 16 - dual task incongruent switch incorrect

eventlist = EEG.event;
eventlist(end+1).type = 'S  1';
stimidx = find(strcmpi({eventlist.type},'S  1'));
handcond = subcondlist(1,2);
EEGold = EEG;

%%
for stimi = 1:length(stimidx)-1

    %% left-hand response
    if any(strcmpi({eventlist(stimidx(stimi):stimidx(stimi+1)-1).type},'L  1'))
        % for handcond 1
        if handcond == 1
            % for red cue
            if ismember(eventlist(stimidx(stimi)).cue,[1,3])
                if ismember(eventlist(stimidx(stimi)).stim,[1,3])
                    eventlist(stimidx(stimi)).corr = 1;
                    respidx = stimidx(stimi) + find(strcmpi({eventlist(stimidx(stimi):stimidx(stimi+1)-1).type},'L  1')) -1 ;
                    respidx = respidx(1);
                    eventlist(stimidx(stimi)).rt = (eventlist(respidx).latency - eventlist(stimidx(stimi)).latency)*1000/EEG.srate;
                else eventlist(stimidx(stimi)).corr = 0;
                    respidx = stimidx(stimi) + find(strcmpi({eventlist(stimidx(stimi):stimidx(stimi+1)-1).type},'L  1')) -1 ;
                    respidx = respidx(1);
                    eventlist(stimidx(stimi)).rt = (eventlist(respidx).latency - eventlist(stimidx(stimi)).latency)*1000/EEG.srate;                end % if ismember
            % for green cue
            elseif ismember(eventlist(stimidx(stimi)).cue,[2,4])
                if ismember(eventlist(stimidx(stimi)).stim,[2,4])
                    eventlist(stimidx(stimi)).corr = 1;
                    respidx = stimidx(stimi) + find(strcmpi({eventlist(stimidx(stimi):stimidx(stimi+1)-1).type},'L  1')) -1 ;
                    respidx = respidx(1);
                    eventlist(stimidx(stimi)).rt = (eventlist(respidx).latency - eventlist(stimidx(stimi)).latency)*1000/EEG.srate;
                else eventlist(stimidx(stimi)).corr = 0;
                    respidx = stimidx(stimi) + find(strcmpi({eventlist(stimidx(stimi):stimidx(stimi+1)-1).type},'L  1')) -1 ;
                    respidx = respidx(1);
                    eventlist(stimidx(stimi)).rt = (eventlist(respidx).latency - eventlist(stimidx(stimi)).latency)*1000/EEG.srate;
                end % if ismember
            end % end cue
        % for handcond 2
        elseif handcond == 2
            % for red cue
            if ismember(eventlist(stimidx(stimi)).cue,[1,3])
                if ismember(eventlist(stimidx(stimi)).stim,[2,4])
                    eventlist(stimidx(stimi)).corr = 1;
                    respidx = stimidx(stimi) + find(strcmpi({eventlist(stimidx(stimi):stimidx(stimi+1)-1).type},'L  1')) -1 ;
                    respidx = respidx(1);
                    eventlist(stimidx(stimi)).rt = (eventlist(respidx).latency - eventlist(stimidx(stimi)).latency)*1000/EEG.srate;
                else eventlist(stimidx(stimi)).corr = 0;
                    respidx = stimidx(stimi) + find(strcmpi({eventlist(stimidx(stimi):stimidx(stimi+1)-1).type},'L  1')) -1 ;
                    respidx = respidx(1);
                    eventlist(stimidx(stimi)).rt = (eventlist(respidx).latency - eventlist(stimidx(stimi)).latency)*1000/EEG.srate;
                end % if ismember
            % for green cue
            elseif ismember(eventlist(stimidx(stimi)).cue,[2,4])
                if ismember(eventlist(stimidx(stimi)).stim,[1,3])
                    eventlist(stimidx(stimi)).corr = 1;
                    respidx = stimidx(stimi) + find(strcmpi({eventlist(stimidx(stimi):stimidx(stimi+1)-1).type},'L  1')) -1 ;
                    respidx = respidx(1);
                    eventlist(stimidx(stimi)).rt = (eventlist(respidx).latency - eventlist(stimidx(stimi)).latency)*1000/EEG.srate;
                else eventlist(stimidx(stimi)).corr = 0;
                    respidx = stimidx(stimi) + find(strcmpi({eventlist(stimidx(stimi):stimidx(stimi+1)-1).type},'L  1')) -1 ;
                    respidx = respidx(1);
                    eventlist(stimidx(stimi)).rt = (eventlist(respidx).latency - eventlist(stimidx(stimi)).latency)*1000/EEG.srate;
                end % if ismember
            end % end cue
        end % if handcond

    %% right hand response
    elseif any(strcmpi({eventlist(stimidx(stimi):stimidx(stimi+1)-1).type},'R  1'))
        % for handcond 1
        if handcond == 1
            % for green cue
            if ismember(eventlist(stimidx(stimi)).cue,[1,3])
                if ismember(eventlist(stimidx(stimi)).stim,[2,4])
                    eventlist(stimidx(stimi)).corr = 1;
                    respidx = stimidx(stimi) + find(strcmpi({eventlist(stimidx(stimi):stimidx(stimi+1)-1).type},'R  1')) -1 ;
                    respidx = respidx(1);
                    eventlist(stimidx(stimi)).rt = (eventlist(respidx).latency - eventlist(stimidx(stimi)).latency)*1000/EEG.srate;
                else eventlist(stimidx(stimi)).corr = 0;
                    respidx = stimidx(stimi) + find(strcmpi({eventlist(stimidx(stimi):stimidx(stimi+1)-1).type},'R  1')) -1 ;
                    respidx = respidx(1);
                    eventlist(stimidx(stimi)).rt = (eventlist(respidx).latency - eventlist(stimidx(stimi)).latency)*1000/EEG.srate;
                end % if ismember
            % for red cue
            elseif ismember(eventlist(stimidx(stimi)).cue,[2,4])
                if ismember(eventlist(stimidx(stimi)).stim,[1,3])
                    eventlist(stimidx(stimi)).corr = 1
                    respidx = stimidx(stimi) + find(strcmpi({eventlist(stimidx(stimi):stimidx(stimi+1)-1).type},'R  1')) -1 ;
                    respidx = respidx(1);
                    eventlist(stimidx(stimi)).rt = (eventlist(respidx).latency - eventlist(stimidx(stimi)).latency)*1000/EEG.srate;
                else eventlist(stimidx(stimi)).corr = 0
                    respidx = stimidx(stimi) + find(strcmpi({eventlist(stimidx(stimi):stimidx(stimi+1)-1).type},'R  1')) -1 ;
                    respidx = respidx(1);
                    eventlist(stimidx(stimi)).rt = (eventlist(respidx).latency - eventlist(stimidx(stimi)).latency)*1000/EEG.srate;
                end % if ismember
            end % end cue
        % for handcond 2
        elseif handcond == 2
            % for green cue
            if ismember(eventlist(stimidx(stimi)).cue,[1,3])
                if ismember(eventlist(stimidx(stimi)).stim,[1,3])
                    eventlist(stimidx(stimi)).corr = 1
                    respidx = stimidx(stimi) + find(strcmpi({eventlist(stimidx(stimi):stimidx(stimi+1)-1).type},'R  1')) -1 ;
                    respidx = respidx(1);
                    eventlist(stimidx(stimi)).rt = (eventlist(respidx).latency - eventlist(stimidx(stimi)).latency)*1000/EEG.srate;
                else eventlist(stimidx(stimi)).corr = 0
                    respidx = stimidx(stimi) + find(strcmpi({eventlist(stimidx(stimi):stimidx(stimi+1)-1).type},'R  1')) -1 ;
                    respidx = respidx(1);
                    eventlist(stimidx(stimi)).rt = (eventlist(respidx).latency - eventlist(stimidx(stimi)).latency)*1000/EEG.srate;
                end % if ismember
            % for red cue
            elseif ismember(eventlist(stimidx(stimi)).cue,[2,4])
                if ismember(eventlist(stimidx(stimi)).stim,[2,4])
                    eventlist(stimidx(stimi)).corr = 1
                    respidx = stimidx(stimi) + find(strcmpi({eventlist(stimidx(stimi):stimidx(stimi+1)-1).type},'R  1')) -1 ;
                    respidx = respidx(1);
                    eventlist(stimidx(stimi)).rt = (eventlist(respidx).latency - eventlist(stimidx(stimi)).latency)*1000/EEG.srate;
                else eventlist(stimidx(stimi)).corr = 0
                    respidx = stimidx(stimi) + find(strcmpi({eventlist(stimidx(stimi):stimidx(stimi+1)-1).type},'R  1')) -1 ;
                    respidx = respidx(1);
                    eventlist(stimidx(stimi)).rt = (eventlist(respidx).latency - eventlist(stimidx(stimi)).latency)*1000/EEG.srate;
                end % if ismember
            end % end cue
        end % if handcond
    else eventlist(stimidx(stimi)).corr = 98;
        reac = 1;
    end % if strcmpi

end % for stimi

eventlist(stimidx(end)) = [];
stimidx = find(strcmpi({eventlist.type},'S  1'));

% set early and late answers to incorrect
for stimi = 1:length(stimidx)
    if eventlist(stimidx(stimi)).rt < 250
        eventlist(stimidx(stimi)).corr = 99;
    elseif eventlist(stimidx(stimi)).rt > 1000
        eventlist(stimidx(stimi)).corr = 99;
    end % if eventlist
end % for stimi

%% calculate mean response times and accuracies
% reduce eventlist
eventlist_redux_idx = find(strcmpi({eventlist.type},'S  1'));
eventlist_r = eventlist(eventlist_redux_idx);

%% auditory
% auditory, single task, congruent
st_a_s_c = find([eventlist_r.move] == 1 & [eventlist_r.modality] == 2 & [eventlist_r.mixed] == 1 & ismember([eventlist_r.stim],[1,2]) & [eventlist_r.corr] == 1);
wa_a_s_c = find([eventlist_r.move] == 2 & [eventlist_r.modality] == 2 & [eventlist_r.mixed] == 1 & ismember([eventlist_r.stim],[1,2]) & [eventlist_r.corr] == 1);
pe_a_s_c = find([eventlist_r.move] == 3 & [eventlist_r.modality] == 2 & [eventlist_r.mixed] == 1 & ismember([eventlist_r.stim],[1,2]) & [eventlist_r.corr] == 1);

% auditory, single task, incongruent
st_a_s_i = find([eventlist_r.move] == 1 & [eventlist_r.modality] == 2 & [eventlist_r.mixed] == 1 & ismember([eventlist_r.stim],[3,4]) & [eventlist_r.corr] == 1);
wa_a_s_i = find([eventlist_r.move] == 2 & [eventlist_r.modality] == 2 & [eventlist_r.mixed] == 1 & ismember([eventlist_r.stim],[3,4]) & [eventlist_r.corr] == 1);
pe_a_s_i = find([eventlist_r.move] == 3 & [eventlist_r.modality] == 2 & [eventlist_r.mixed] == 1 & ismember([eventlist_r.stim],[3,4]) & [eventlist_r.corr] == 1);

% auditory, dual task, congruent
st_a_d_c = find([eventlist_r.move] == 1 & [eventlist_r.modality] == 2 & [eventlist_r.mixed] == 2 & ismember([eventlist_r.stim],[1,2]) & [eventlist_r.corr] == 1);
wa_a_d_c = find([eventlist_r.move] == 2 & [eventlist_r.modality] == 2 & [eventlist_r.mixed] == 2 & ismember([eventlist_r.stim],[1,2]) & [eventlist_r.corr] == 1);
pe_a_d_c = find([eventlist_r.move] == 3 & [eventlist_r.modality] == 2 & [eventlist_r.mixed] == 2 & ismember([eventlist_r.stim],[1,2]) & [eventlist_r.corr] == 1);

% auditory, dual task, incongruent
st_a_d_i = find([eventlist_r.move] == 1 & [eventlist_r.modality] == 2 & [eventlist_r.mixed] == 2 & ismember([eventlist_r.stim],[3,4]) & [eventlist_r.corr] == 1);
wa_a_d_i = find([eventlist_r.move] == 2 & [eventlist_r.modality] == 2 & [eventlist_r.mixed] == 2 & ismember([eventlist_r.stim],[3,4]) & [eventlist_r.corr] == 1);
pe_a_d_i = find([eventlist_r.move] == 3 & [eventlist_r.modality] == 2 & [eventlist_r.mixed] == 2 & ismember([eventlist_r.stim],[3,4]) & [eventlist_r.corr] == 1);

%% visual
% visual, single task, congruent, repeat
st_v_s_c = find([eventlist_r.move] == 1 & [eventlist_r.modality] == 1 & [eventlist_r.mixed] == 1 & ismember([eventlist_r.stim],[1,2]) & [eventlist_r.corr] == 1);
wa_v_s_c = find([eventlist_r.move] == 2 & [eventlist_r.modality] == 1 & [eventlist_r.mixed] == 1 & ismember([eventlist_r.stim],[1,2]) & [eventlist_r.corr] == 1);
pe_v_s_c = find([eventlist_r.move] == 3 & [eventlist_r.modality] == 1 & [eventlist_r.mixed] == 1 & ismember([eventlist_r.stim],[1,2]) & [eventlist_r.corr] == 1);

% auditory, single task, incongruent
st_v_s_i = find([eventlist_r.move] == 1 & [eventlist_r.modality] == 1 & [eventlist_r.mixed] == 1 & ismember([eventlist_r.stim],[3,4]) & [eventlist_r.corr] == 1);
wa_v_s_i = find([eventlist_r.move] == 2 & [eventlist_r.modality] == 1 & [eventlist_r.mixed] == 1 & ismember([eventlist_r.stim],[3,4]) & [eventlist_r.corr] == 1);
pe_v_s_i = find([eventlist_r.move] == 3 & [eventlist_r.modality] == 1 & [eventlist_r.mixed] == 1 & ismember([eventlist_r.stim],[3,4]) & [eventlist_r.corr] == 1);

% auditory, dual task, congruent
st_v_d_c = find([eventlist_r.move] == 1 & [eventlist_r.modality] == 1 & [eventlist_r.mixed] == 2 & ismember([eventlist_r.stim],[1,2]) & [eventlist_r.corr] == 1);
wa_v_d_c = find([eventlist_r.move] == 2 & [eventlist_r.modality] == 1 & [eventlist_r.mixed] == 2 & ismember([eventlist_r.stim],[1,2]) & [eventlist_r.corr] == 1);
pe_v_d_c = find([eventlist_r.move] == 3 & [eventlist_r.modality] == 1 & [eventlist_r.mixed] == 2 & ismember([eventlist_r.stim],[1,2]) & [eventlist_r.corr] == 1);

% auditory, dual task, incongruent
st_v_d_i = find([eventlist_r.move] == 1 & [eventlist_r.modality] == 1 & [eventlist_r.mixed] == 2 & ismember([eventlist_r.stim],[3,4]) & [eventlist_r.corr] == 1);
wa_v_d_i = find([eventlist_r.move] == 2 & [eventlist_r.modality] == 1 & [eventlist_r.mixed] == 2 & ismember([eventlist_r.stim],[3,4]) & [eventlist_r.corr] == 1);
pe_v_d_i = find([eventlist_r.move] == 3 & [eventlist_r.modality] == 1 & [eventlist_r.mixed] == 2 & ismember([eventlist_r.stim],[3,4]) & [eventlist_r.corr] == 1);

% index all conditions
% find stay and switch trials
for switchi = 2:length(eventlist_r)
    lastCue = eventlist_r(switchi-1).cue;
    currCue = eventlist_r(switchi).cue;
    if lastCue == currCue
        eventlist_r(switchi).switch = 0;
    elseif lastCue ~= currCue
        eventlist_r(switchi).switch = 1;
    end % if lastCue
end % for switchi

% re-assign first event of each block as repeat
eventlist_r(st_a_s_c(1)).switch = 0;
eventlist_r(wa_a_s_c(1)).switch = 0;
eventlist_r(pe_a_s_c(1)).switch = 0;
eventlist_r(st_a_s_i(1)).switch = 0;
eventlist_r(wa_a_s_i(1)).switch = 0;
eventlist_r(pe_a_s_i(1)).switch = 0;
eventlist_r(st_a_d_c(1)).switch = 0;
eventlist_r(wa_a_d_c(1)).switch = 0;
eventlist_r(pe_a_d_c(1)).switch = 0;
eventlist_r(st_a_d_i(1)).switch = 0;
eventlist_r(wa_a_d_i(1)).switch = 0;
eventlist_r(pe_a_d_i(1)).switch = 0;

eventlist_r(st_v_s_c(1)).switch = 0;
eventlist_r(wa_v_s_c(1)).switch = 0;
eventlist_r(pe_v_s_c(1)).switch = 0;
eventlist_r(st_v_s_i(1)).switch = 0;
eventlist_r(wa_v_s_i(1)).switch = 0;
eventlist_r(pe_v_s_i(1)).switch = 0;
eventlist_r(st_v_d_c(1)).switch = 0;
eventlist_r(wa_v_d_c(1)).switch = 0;
eventlist_r(pe_v_d_c(1)).switch = 0;
eventlist_r(st_v_d_i(1)).switch = 0;
eventlist_r(wa_v_d_i(1)).switch = 0;
eventlist_r(pe_v_d_i(1)).switch = 0;

eventlist_r(1).switch = 0;
%%
% auditory, dual task, congruent, repeat
st_a_d_c_r = find([eventlist_r.move] == 1 & [eventlist_r.modality] == 2 & [eventlist_r.mixed] == 2 & ismember([eventlist_r.stim],[1,2]) & [eventlist_r.corr] == 1 & [eventlist_r.switch] == 0);
wa_a_d_c_r = find([eventlist_r.move] == 2 & [eventlist_r.modality] == 2 & [eventlist_r.mixed] == 2 & ismember([eventlist_r.stim],[1,2]) & [eventlist_r.corr] == 1 & [eventlist_r.switch] == 0);
pe_a_d_c_r = find([eventlist_r.move] == 3 & [eventlist_r.modality] == 2 & [eventlist_r.mixed] == 2 & ismember([eventlist_r.stim],[1,2]) & [eventlist_r.corr] == 1 & [eventlist_r.switch] == 0);

% auditory, dual task, congruent, switch
st_a_d_c_s = find([eventlist_r.move] == 1 & [eventlist_r.modality] == 2 & [eventlist_r.mixed] == 2 & ismember([eventlist_r.stim],[1,2]) & [eventlist_r.corr] == 1 & [eventlist_r.switch] == 1);
wa_a_d_c_s = find([eventlist_r.move] == 2 & [eventlist_r.modality] == 2 & [eventlist_r.mixed] == 2 & ismember([eventlist_r.stim],[1,2]) & [eventlist_r.corr] == 1 & [eventlist_r.switch] == 1);
pe_a_d_c_s = find([eventlist_r.move] == 3 & [eventlist_r.modality] == 2 & [eventlist_r.mixed] == 2 & ismember([eventlist_r.stim],[1,2]) & [eventlist_r.corr] == 1 & [eventlist_r.switch] == 1);

% auditory, dual task, incongruent, repeat
st_a_d_i_r = find([eventlist_r.move] == 1 & [eventlist_r.modality] == 2 & [eventlist_r.mixed] == 2 & ismember([eventlist_r.stim],[3,4]) & [eventlist_r.corr] == 1 & [eventlist_r.switch] == 0);
wa_a_d_i_r = find([eventlist_r.move] == 2 & [eventlist_r.modality] == 2 & [eventlist_r.mixed] == 2 & ismember([eventlist_r.stim],[3,4]) & [eventlist_r.corr] == 1 & [eventlist_r.switch] == 0);
pe_a_d_i_r = find([eventlist_r.move] == 3 & [eventlist_r.modality] == 2 & [eventlist_r.mixed] == 2 & ismember([eventlist_r.stim],[3,4]) & [eventlist_r.corr] == 1 & [eventlist_r.switch] == 0);

% auditory, dual task, incongruent, switch
st_a_d_i_s = find([eventlist_r.move] == 1 & [eventlist_r.modality] == 2 & [eventlist_r.mixed] == 2 & ismember([eventlist_r.stim],[3,4]) & [eventlist_r.corr] == 1 & [eventlist_r.switch] == 1);
wa_a_d_i_s = find([eventlist_r.move] == 2 & [eventlist_r.modality] == 2 & [eventlist_r.mixed] == 2 & ismember([eventlist_r.stim],[3,4]) & [eventlist_r.corr] == 1 & [eventlist_r.switch] == 1);
pe_a_d_i_s = find([eventlist_r.move] == 3 & [eventlist_r.modality] == 2 & [eventlist_r.mixed] == 2 & ismember([eventlist_r.stim],[3,4]) & [eventlist_r.corr] == 1 & [eventlist_r.switch] == 1);

%% incorrect auditory
% auditory, single task, congruent
st_a_s_c_incor = find([eventlist_r.move] == 1 & [eventlist_r.modality] == 2 & [eventlist_r.mixed] == 1 & ismember([eventlist_r.stim],[1,2]) & [eventlist_r.corr] == 0);
wa_a_s_c_incor = find([eventlist_r.move] == 2 & [eventlist_r.modality] == 2 & [eventlist_r.mixed] == 1 & ismember([eventlist_r.stim],[1,2]) & [eventlist_r.corr] == 0);
pe_a_s_c_incor = find([eventlist_r.move] == 3 & [eventlist_r.modality] == 2 & [eventlist_r.mixed] == 1 & ismember([eventlist_r.stim],[1,2]) & [eventlist_r.corr] == 0);

% auditory, single task, incongruent
st_a_s_i_incor = find([eventlist_r.move] == 1 & [eventlist_r.modality] == 2 & [eventlist_r.mixed] == 1 & ismember([eventlist_r.stim],[3,4]) & [eventlist_r.corr] == 0);
wa_a_s_i_incor = find([eventlist_r.move] == 2 & [eventlist_r.modality] == 2 & [eventlist_r.mixed] == 1 & ismember([eventlist_r.stim],[3,4]) & [eventlist_r.corr] == 0);
pe_a_s_i_incor = find([eventlist_r.move] == 3 & [eventlist_r.modality] == 2 & [eventlist_r.mixed] == 1 & ismember([eventlist_r.stim],[3,4]) & [eventlist_r.corr] == 0);

% auditory, dual task, congruent
st_a_d_c_incor = find([eventlist_r.move] == 1 & [eventlist_r.modality] == 2 & [eventlist_r.mixed] == 2 & ismember([eventlist_r.stim],[1,2]) & [eventlist_r.corr] == 0);
wa_a_d_c_incor = find([eventlist_r.move] == 2 & [eventlist_r.modality] == 2 & [eventlist_r.mixed] == 2 & ismember([eventlist_r.stim],[1,2]) & [eventlist_r.corr] == 0);
pe_a_d_c_incor = find([eventlist_r.move] == 3 & [eventlist_r.modality] == 2 & [eventlist_r.mixed] == 2 & ismember([eventlist_r.stim],[1,2]) & [eventlist_r.corr] == 0);

% auditory, dual task, incongruent
st_a_d_i_incor = find([eventlist_r.move] == 1 & [eventlist_r.modality] == 2 & [eventlist_r.mixed] == 2 & ismember([eventlist_r.stim],[3,4]) & [eventlist_r.corr] == 0);
wa_a_d_i_incor = find([eventlist_r.move] == 2 & [eventlist_r.modality] == 2 & [eventlist_r.mixed] == 2 & ismember([eventlist_r.stim],[3,4]) & [eventlist_r.corr] == 0);
pe_a_d_i_incor = find([eventlist_r.move] == 3 & [eventlist_r.modality] == 2 & [eventlist_r.mixed] == 2 & ismember([eventlist_r.stim],[3,4]) & [eventlist_r.corr] == 0);

% auditory, dual task, congruent, repeat
st_a_d_c_r_incor = find([eventlist_r.move] == 1 & [eventlist_r.modality] == 2 & [eventlist_r.mixed] == 2 & ismember([eventlist_r.stim],[1,2]) & [eventlist_r.corr] == 0 & [eventlist_r.switch] == 0);
wa_a_d_c_r_incor = find([eventlist_r.move] == 2 & [eventlist_r.modality] == 2 & [eventlist_r.mixed] == 2 & ismember([eventlist_r.stim],[1,2]) & [eventlist_r.corr] == 0 & [eventlist_r.switch] == 0);
pe_a_d_c_r_incor = find([eventlist_r.move] == 3 & [eventlist_r.modality] == 2 & [eventlist_r.mixed] == 2 & ismember([eventlist_r.stim],[1,2]) & [eventlist_r.corr] == 0 & [eventlist_r.switch] == 0);

% auditory, dual task, congruent, switch
st_a_d_c_s_incor = find([eventlist_r.move] == 1 & [eventlist_r.modality] == 2 & [eventlist_r.mixed] == 2 & ismember([eventlist_r.stim],[1,2]) & [eventlist_r.corr] == 0 & [eventlist_r.switch] == 1);
wa_a_d_c_s_incor = find([eventlist_r.move] == 2 & [eventlist_r.modality] == 2 & [eventlist_r.mixed] == 2 & ismember([eventlist_r.stim],[1,2]) & [eventlist_r.corr] == 0 & [eventlist_r.switch] == 1);
pe_a_d_c_s_incor = find([eventlist_r.move] == 3 & [eventlist_r.modality] == 2 & [eventlist_r.mixed] == 2 & ismember([eventlist_r.stim],[1,2]) & [eventlist_r.corr] == 0 & [eventlist_r.switch] == 1);

% auditory, dual task, incongruent, repeat
st_a_d_i_r_incor = find([eventlist_r.move] == 1 & [eventlist_r.modality] == 2 & [eventlist_r.mixed] == 2 & ismember([eventlist_r.stim],[3,4]) & [eventlist_r.corr] == 0 & [eventlist_r.switch] == 0);
wa_a_d_i_r_incor = find([eventlist_r.move] == 2 & [eventlist_r.modality] == 2 & [eventlist_r.mixed] == 2 & ismember([eventlist_r.stim],[3,4]) & [eventlist_r.corr] == 0 & [eventlist_r.switch] == 0);
pe_a_d_i_r_incor = find([eventlist_r.move] == 3 & [eventlist_r.modality] == 2 & [eventlist_r.mixed] == 2 & ismember([eventlist_r.stim],[3,4]) & [eventlist_r.corr] == 0 & [eventlist_r.switch] == 0);

% auditory, dual task, incongruent, switch
st_a_d_i_s_incor = find([eventlist_r.move] == 1 & [eventlist_r.modality] == 2 & [eventlist_r.mixed] == 2 & ismember([eventlist_r.stim],[3,4]) & [eventlist_r.corr] == 0 & [eventlist_r.switch] == 1);
wa_a_d_i_s_incor = find([eventlist_r.move] == 2 & [eventlist_r.modality] == 2 & [eventlist_r.mixed] == 2 & ismember([eventlist_r.stim],[3,4]) & [eventlist_r.corr] == 0 & [eventlist_r.switch] == 1);
pe_a_d_i_s_incor = find([eventlist_r.move] == 3 & [eventlist_r.modality] == 2 & [eventlist_r.mixed] == 2 & ismember([eventlist_r.stim],[3,4]) & [eventlist_r.corr] == 0 & [eventlist_r.switch] == 1);

%%
% auditory, dual task, congruent, repeat
st_v_d_c_r = find([eventlist_r.move] == 1 & [eventlist_r.modality] == 1 & [eventlist_r.mixed] == 2 & ismember([eventlist_r.stim],[1,2]) & [eventlist_r.corr] == 1 & [eventlist_r.switch] == 0);
wa_v_d_c_r = find([eventlist_r.move] == 2 & [eventlist_r.modality] == 1 & [eventlist_r.mixed] == 2 & ismember([eventlist_r.stim],[1,2]) & [eventlist_r.corr] == 1 & [eventlist_r.switch] == 0);
pe_v_d_c_r = find([eventlist_r.move] == 3 & [eventlist_r.modality] == 1 & [eventlist_r.mixed] == 2 & ismember([eventlist_r.stim],[1,2]) & [eventlist_r.corr] == 1 & [eventlist_r.switch] == 0);

% auditory, dual task, congruent, switch
st_v_d_c_s = find([eventlist_r.move] == 1 & [eventlist_r.modality] == 1 & [eventlist_r.mixed] == 2 & ismember([eventlist_r.stim],[1,2]) & [eventlist_r.corr] == 1 & [eventlist_r.switch] == 1);
wa_v_d_c_s = find([eventlist_r.move] == 2 & [eventlist_r.modality] == 1 & [eventlist_r.mixed] == 2 & ismember([eventlist_r.stim],[1,2]) & [eventlist_r.corr] == 1 & [eventlist_r.switch] == 1);
pe_v_d_c_s = find([eventlist_r.move] == 3 & [eventlist_r.modality] == 1 & [eventlist_r.mixed] == 2 & ismember([eventlist_r.stim],[1,2]) & [eventlist_r.corr] == 1 & [eventlist_r.switch] == 1);

% auditory, dual task, incongruent, repeat
st_v_d_i_r = find([eventlist_r.move] == 1 & [eventlist_r.modality] == 1 & [eventlist_r.mixed] == 2 & ismember([eventlist_r.stim],[3,4]) & [eventlist_r.corr] == 1 & [eventlist_r.switch] == 0);
wa_v_d_i_r = find([eventlist_r.move] == 2 & [eventlist_r.modality] == 1 & [eventlist_r.mixed] == 2 & ismember([eventlist_r.stim],[3,4]) & [eventlist_r.corr] == 1 & [eventlist_r.switch] == 0);
pe_v_d_i_r = find([eventlist_r.move] == 3 & [eventlist_r.modality] == 1 & [eventlist_r.mixed] == 2 & ismember([eventlist_r.stim],[3,4]) & [eventlist_r.corr] == 1 & [eventlist_r.switch] == 0);

% auditory, dual task, incongruent, switch
st_v_d_i_s = find([eventlist_r.move] == 1 & [eventlist_r.modality] == 1 & [eventlist_r.mixed] == 2 & ismember([eventlist_r.stim],[3,4]) & [eventlist_r.corr] == 1 & [eventlist_r.switch] == 1);
wa_v_d_i_s = find([eventlist_r.move] == 2 & [eventlist_r.modality] == 1 & [eventlist_r.mixed] == 2 & ismember([eventlist_r.stim],[3,4]) & [eventlist_r.corr] == 1 & [eventlist_r.switch] == 1);
pe_v_d_i_s = find([eventlist_r.move] == 3 & [eventlist_r.modality] == 1 & [eventlist_r.mixed] == 2 & ismember([eventlist_r.stim],[3,4]) & [eventlist_r.corr] == 1 & [eventlist_r.switch] == 1);

%% incorrect visual
% visual, single task, congruent, repeat
st_v_s_c_incor = find([eventlist_r.move] == 1 & [eventlist_r.modality] == 1 & [eventlist_r.mixed] == 1 & ismember([eventlist_r.stim],[1,2]) & [eventlist_r.corr] == 0);
wa_v_s_c_incor = find([eventlist_r.move] == 2 & [eventlist_r.modality] == 1 & [eventlist_r.mixed] == 1 & ismember([eventlist_r.stim],[1,2]) & [eventlist_r.corr] == 0);
pe_v_s_c_incor = find([eventlist_r.move] == 3 & [eventlist_r.modality] == 1 & [eventlist_r.mixed] == 1 & ismember([eventlist_r.stim],[1,2]) & [eventlist_r.corr] == 0);

% auditory, single task, incongruent
st_v_s_i_incor = find([eventlist_r.move] == 1 & [eventlist_r.modality] == 1 & [eventlist_r.mixed] == 1 & ismember([eventlist_r.stim],[3,4]) & [eventlist_r.corr] == 0);
wa_v_s_i_incor = find([eventlist_r.move] == 2 & [eventlist_r.modality] == 1 & [eventlist_r.mixed] == 1 & ismember([eventlist_r.stim],[3,4]) & [eventlist_r.corr] == 0);
pe_v_s_i_incor = find([eventlist_r.move] == 3 & [eventlist_r.modality] == 1 & [eventlist_r.mixed] == 1 & ismember([eventlist_r.stim],[3,4]) & [eventlist_r.corr] == 0);

% auditory, dual task, congruent
st_v_d_c_incor = find([eventlist_r.move] == 1 & [eventlist_r.modality] == 1 & [eventlist_r.mixed] == 2 & ismember([eventlist_r.stim],[1,2]) & [eventlist_r.corr] == 0);
wa_v_d_c_incor = find([eventlist_r.move] == 2 & [eventlist_r.modality] == 1 & [eventlist_r.mixed] == 2 & ismember([eventlist_r.stim],[1,2]) & [eventlist_r.corr] == 0);
pe_v_d_c_incor = find([eventlist_r.move] == 3 & [eventlist_r.modality] == 1 & [eventlist_r.mixed] == 2 & ismember([eventlist_r.stim],[1,2]) & [eventlist_r.corr] == 0);

% auditory, dual task, incongruent
st_v_d_i_incor = find([eventlist_r.move] == 1 & [eventlist_r.modality] == 1 & [eventlist_r.mixed] == 2 & ismember([eventlist_r.stim],[3,4]) & [eventlist_r.corr] == 0);
wa_v_d_i_incor = find([eventlist_r.move] == 2 & [eventlist_r.modality] == 1 & [eventlist_r.mixed] == 2 & ismember([eventlist_r.stim],[3,4]) & [eventlist_r.corr] == 0);
pe_v_d_i_incor = find([eventlist_r.move] == 3 & [eventlist_r.modality] == 1 & [eventlist_r.mixed] == 2 & ismember([eventlist_r.stim],[3,4]) & [eventlist_r.corr] == 0);

% auditory, dual task, congruent, repeat
st_v_d_c_r_incor = find([eventlist_r.move] == 1 & [eventlist_r.modality] == 1 & [eventlist_r.mixed] == 2 & ismember([eventlist_r.stim],[1,2]) & [eventlist_r.corr] == 0 & [eventlist_r.switch] == 0);
wa_v_d_c_r_incor = find([eventlist_r.move] == 2 & [eventlist_r.modality] == 1 & [eventlist_r.mixed] == 2 & ismember([eventlist_r.stim],[1,2]) & [eventlist_r.corr] == 0 & [eventlist_r.switch] == 0);
pe_v_d_c_r_incor = find([eventlist_r.move] == 3 & [eventlist_r.modality] == 1 & [eventlist_r.mixed] == 2 & ismember([eventlist_r.stim],[1,2]) & [eventlist_r.corr] == 0 & [eventlist_r.switch] == 0);

% auditory, dual task, congruent, switch
st_v_d_c_s_incor = find([eventlist_r.move] == 1 & [eventlist_r.modality] == 1 & [eventlist_r.mixed] == 2 & ismember([eventlist_r.stim],[1,2]) & [eventlist_r.corr] == 0 & [eventlist_r.switch] == 1);
wa_v_d_c_s_incor = find([eventlist_r.move] == 2 & [eventlist_r.modality] == 1 & [eventlist_r.mixed] == 2 & ismember([eventlist_r.stim],[1,2]) & [eventlist_r.corr] == 0 & [eventlist_r.switch] == 1);
pe_v_d_c_s_incor = find([eventlist_r.move] == 3 & [eventlist_r.modality] == 1 & [eventlist_r.mixed] == 2 & ismember([eventlist_r.stim],[1,2]) & [eventlist_r.corr] == 0 & [eventlist_r.switch] == 1);

% auditory, dual task, incongruent, repeat
st_v_d_i_r_incor = find([eventlist_r.move] == 1 & [eventlist_r.modality] == 1 & [eventlist_r.mixed] == 2 & ismember([eventlist_r.stim],[3,4]) & [eventlist_r.corr] == 0 & [eventlist_r.switch] == 0);
wa_v_d_i_r_incor = find([eventlist_r.move] == 2 & [eventlist_r.modality] == 1 & [eventlist_r.mixed] == 2 & ismember([eventlist_r.stim],[3,4]) & [eventlist_r.corr] == 0 & [eventlist_r.switch] == 0);
pe_v_d_i_r_incor = find([eventlist_r.move] == 3 & [eventlist_r.modality] == 1 & [eventlist_r.mixed] == 2 & ismember([eventlist_r.stim],[3,4]) & [eventlist_r.corr] == 0 & [eventlist_r.switch] == 0);

% auditory, dual task, incongruent, switch
st_v_d_i_s_incor = find([eventlist_r.move] == 1 & [eventlist_r.modality] == 1 & [eventlist_r.mixed] == 2 & ismember([eventlist_r.stim],[3,4]) & [eventlist_r.corr] == 0 & [eventlist_r.switch] == 1);
wa_v_d_i_s_incor = find([eventlist_r.move] == 2 & [eventlist_r.modality] == 1 & [eventlist_r.mixed] == 2 & ismember([eventlist_r.stim],[3,4]) & [eventlist_r.corr] == 0 & [eventlist_r.switch] == 1);
pe_v_d_i_s_incor = find([eventlist_r.move] == 3 & [eventlist_r.modality] == 1 & [eventlist_r.mixed] == 2 & ismember([eventlist_r.stim],[3,4]) & [eventlist_r.corr] == 0 & [eventlist_r.switch] == 1);

%% calculate
reacmat_a = []; reacmat_v = []; accumat_a = []; accumat_v = [];

% description for cols (rows always standt,walk,pert):
%  1 - single task congruent
%  2 - single task incongruent
%  3 - dual task congruent
%  4 - dual task incongruent
%  5 - dual task congruent repeat
%  6 - dual task incongruent repeat
%  7 - dual task congruent switch
%  8 - dual task incongruent switch
%  9 - single task congruent incorrect
% 10 - single task incongruent incorrect
% 11 - dual task congruent incorrect
% 12 - dual task incongruent incorrect
% 13 - dual task congruent repeat incorrect
% 14 - dual task incongruent repeat incorrect
% 15 - dual task congruent switch incorrect
% 16 - dual task incongruent switch incorrect

reacmat_a = [mean([eventlist_r(st_a_s_c).rt],'omitnan'), mean([eventlist_r(wa_a_s_c).rt],'omitnan'), mean([eventlist_r(pe_a_s_c).rt],'omitnan');...
             mean([eventlist_r(st_a_s_i).rt],'omitnan'), mean([eventlist_r(wa_a_s_i).rt],'omitnan'), mean([eventlist_r(pe_a_s_i).rt],'omitnan');...
             mean([eventlist_r(st_a_d_c).rt],'omitnan'), mean([eventlist_r(wa_a_d_c).rt],'omitnan'), mean([eventlist_r(pe_a_d_c).rt],'omitnan');...
             mean([eventlist_r(st_a_d_i).rt],'omitnan'), mean([eventlist_r(wa_a_d_i).rt],'omitnan'), mean([eventlist_r(pe_a_d_i).rt],'omitnan');...
             mean([eventlist_r(st_a_d_c_r).rt],'omitnan'), mean([eventlist_r(wa_a_d_c_r).rt],'omitnan'), mean([eventlist_r(pe_a_d_c_r).rt],'omitnan');...
             mean([eventlist_r(st_a_d_i_r).rt],'omitnan'), mean([eventlist_r(wa_a_d_i_r).rt],'omitnan'), mean([eventlist_r(pe_a_d_i_r).rt],'omitnan');...
             mean([eventlist_r(st_a_d_c_s).rt],'omitnan'), mean([eventlist_r(wa_a_d_c_s).rt],'omitnan'), mean([eventlist_r(pe_a_d_c_s).rt],'omitnan');...
             mean([eventlist_r(st_a_d_i_s).rt],'omitnan'), mean([eventlist_r(wa_a_d_i_s).rt],'omitnan'), mean([eventlist_r(pe_a_d_i_s).rt],'omitnan');...
             mean([eventlist_r(st_a_s_c_incor).rt],'omitnan'), mean([eventlist_r(wa_a_s_c_incor).rt],'omitnan'), mean([eventlist_r(pe_a_s_c_incor).rt],'omitnan');...
             mean([eventlist_r(st_a_s_i_incor).rt],'omitnan'), mean([eventlist_r(wa_a_s_i_incor).rt],'omitnan'), mean([eventlist_r(pe_a_s_i_incor).rt],'omitnan');...
             mean([eventlist_r(st_a_d_c_incor).rt],'omitnan'), mean([eventlist_r(wa_a_d_c_incor).rt],'omitnan'), mean([eventlist_r(pe_a_d_c_incor).rt],'omitnan');...
             mean([eventlist_r(st_a_d_i_incor).rt],'omitnan'), mean([eventlist_r(wa_a_d_i_incor).rt],'omitnan'), mean([eventlist_r(pe_a_d_i_incor).rt],'omitnan');...
             mean([eventlist_r(st_a_d_c_r_incor).rt],'omitnan'), mean([eventlist_r(wa_a_d_c_r_incor).rt],'omitnan'), mean([eventlist_r(pe_a_d_c_r_incor).rt],'omitnan');...
             mean([eventlist_r(st_a_d_i_r_incor).rt],'omitnan'), mean([eventlist_r(wa_a_d_i_r_incor).rt],'omitnan'), mean([eventlist_r(pe_a_d_i_r_incor).rt],'omitnan');...
             mean([eventlist_r(st_a_d_c_s_incor).rt],'omitnan'), mean([eventlist_r(wa_a_d_c_s_incor).rt],'omitnan'), mean([eventlist_r(pe_a_d_c_s_incor).rt],'omitnan');...
             mean([eventlist_r(st_a_d_i_s_incor).rt],'omitnan'), mean([eventlist_r(wa_a_d_i_s_incor).rt],'omitnan'), mean([eventlist_r(pe_a_d_i_s_incor).rt],'omitnan')];

reacmat_v = [mean([eventlist_r(st_v_s_c).rt]), mean([eventlist_r(wa_v_s_c).rt]), mean([eventlist_r(pe_v_s_c).rt]);...
             mean([eventlist_r(st_v_s_i).rt]), mean([eventlist_r(wa_v_s_i).rt]), mean([eventlist_r(pe_v_s_i).rt]);...
             mean([eventlist_r(st_v_d_c).rt]), mean([eventlist_r(wa_v_d_c).rt]), mean([eventlist_r(pe_v_d_c).rt]);...
             mean([eventlist_r(st_v_d_i).rt]), mean([eventlist_r(wa_v_d_i).rt]), mean([eventlist_r(pe_v_d_i).rt]);...
             mean([eventlist_r(st_v_d_c_r).rt]), mean([eventlist_r(wa_v_d_c_r).rt]), mean([eventlist_r(pe_v_d_c_r).rt]);...
             mean([eventlist_r(st_v_d_i_r).rt]), mean([eventlist_r(wa_v_d_i_r).rt]), mean([eventlist_r(pe_v_d_i_r).rt]);...
             mean([eventlist_r(st_v_d_c_s).rt]), mean([eventlist_r(wa_v_d_c_s).rt]), mean([eventlist_r(pe_v_d_c_s).rt]);...
             mean([eventlist_r(st_v_d_i_s).rt]), mean([eventlist_r(wa_v_d_i_s).rt]), mean([eventlist_r(pe_v_d_i_s).rt]);...
             mean([eventlist_r(st_v_s_c_incor).rt]), mean([eventlist_r(wa_v_s_c_incor).rt]), mean([eventlist_r(pe_v_s_c_incor).rt]);...
             mean([eventlist_r(st_v_s_i_incor).rt]), mean([eventlist_r(wa_v_s_i_incor).rt]), mean([eventlist_r(pe_v_s_i_incor).rt]);...
             mean([eventlist_r(st_v_d_c_incor).rt]), mean([eventlist_r(wa_v_d_c_incor).rt]), mean([eventlist_r(pe_v_d_c_incor).rt]);...
             mean([eventlist_r(st_v_d_i_incor).rt]), mean([eventlist_r(wa_v_d_i_incor).rt]), mean([eventlist_r(pe_v_d_i_incor).rt]);...
             mean([eventlist_r(st_v_d_c_r_incor).rt]), mean([eventlist_r(wa_v_d_c_r_incor).rt]), mean([eventlist_r(pe_v_d_c_r_incor).rt]);...
             mean([eventlist_r(st_v_d_i_r_incor).rt]), mean([eventlist_r(wa_v_d_i_r_incor).rt]), mean([eventlist_r(pe_v_d_i_r_incor).rt]);...
             mean([eventlist_r(st_v_d_c_s_incor).rt]), mean([eventlist_r(wa_v_d_c_s_incor).rt]), mean([eventlist_r(pe_v_d_c_s_incor).rt]);...
             mean([eventlist_r(st_v_d_i_s_incor).rt]), mean([eventlist_r(wa_v_d_i_s_incor).rt]), mean([eventlist_r(pe_v_d_i_s_incor).rt])];

accumat_a = [length(st_a_s_c)/(length(st_a_s_c)+length(st_a_s_c_incor)), length(wa_a_s_c)/(length(wa_a_s_c)+length(wa_a_s_c_incor)), length(pe_a_s_c)/(length(pe_a_s_c)+length(pe_a_s_c_incor));...
             length(st_a_s_i)/(length(st_a_s_i)+length(st_a_s_i_incor)), length(wa_a_s_i)/(length(wa_a_s_i)+length(wa_a_s_i_incor)), length(pe_a_s_i)/(length(pe_a_s_i)+length(pe_a_s_i_incor));...
             length(st_a_d_c)/(length(st_a_d_c)+length(st_a_d_c_incor)), length(wa_a_d_c)/(length(wa_a_d_c)+length(wa_a_d_c_incor)), length(pe_a_d_c)/(length(pe_a_d_c)+length(pe_a_d_c_incor));...
             length(st_a_d_i)/(length(st_a_d_i)+length(st_a_d_i_incor)), length(wa_a_d_i)/(length(wa_a_d_i)+length(wa_a_d_i_incor)), length(pe_a_d_i)/(length(pe_a_d_i)+length(pe_a_d_i_incor));...
             length(st_a_d_c_r)/(length(st_a_d_c_r)+length(st_a_d_c_r_incor)), length(wa_a_d_c_r)/(length(wa_a_d_c_r)+length(wa_a_d_c_r_incor)), length(pe_a_d_c_r)/(length(pe_a_d_c_r)+length(pe_a_d_c_r_incor));...
             length(st_a_d_i_r)/(length(st_a_d_i_r)+length(st_a_d_i_r_incor)), length(wa_a_d_i_r)/(length(wa_a_d_i_r)+length(wa_a_d_i_r_incor)), length(pe_a_d_i_r)/(length(pe_a_d_i_r)+length(pe_a_d_i_r_incor));...
             length(st_a_d_c_s)/(length(st_a_d_c_s)+length(st_a_d_c_s_incor)), length(wa_a_d_c_s)/(length(wa_a_d_c_s)+length(wa_a_d_c_s_incor)), length(pe_a_d_c_s)/(length(pe_a_d_c_s)+length(pe_a_d_c_s_incor));...
             length(st_a_d_i_s)/(length(st_a_d_i_s)+length(st_a_d_i_s_incor)), length(wa_a_d_i_s)/(length(wa_a_d_i_s)+length(wa_a_d_i_s_incor)), length(pe_a_d_i_s)/(length(pe_a_d_i_s)+length(pe_a_d_i_s_incor))];

accumat_v = [length(st_v_s_c)/(length(st_v_s_c)+length(st_v_s_c_incor)), length(wa_v_s_c)/(length(wa_v_s_c)+length(wa_v_s_c_incor)), length(pe_v_s_c)/(length(pe_v_s_c)+length(pe_v_s_c_incor));...
             length(st_v_s_i)/(length(st_v_s_i)+length(st_v_s_i_incor)), length(wa_v_s_i)/(length(wa_v_s_i)+length(wa_v_s_i_incor)), length(pe_v_s_i)/(length(pe_v_s_i)+length(pe_v_s_i_incor));...
             length(st_v_d_c)/(length(st_v_d_c)+length(st_v_d_c_incor)), length(wa_v_d_c)/(length(wa_v_d_c)+length(wa_v_d_c_incor)), length(pe_v_d_c)/(length(pe_v_d_c)+length(pe_v_d_c_incor));...
             length(st_v_d_i)/(length(st_v_d_i)+length(st_v_d_i_incor)), length(wa_v_d_i)/(length(wa_v_d_i)+length(wa_v_d_i_incor)), length(pe_v_d_i)/(length(pe_v_d_i)+length(pe_v_d_i_incor));...
             length(st_v_d_c_r)/(length(st_v_d_c_r)+length(st_v_d_c_r_incor)), length(wa_v_d_c_r)/(length(wa_v_d_c_r)+length(wa_v_d_c_r_incor)), length(pe_v_d_c_r)/(length(pe_v_d_c_r)+length(pe_v_d_c_r_incor));...
             length(st_v_d_i_r)/(length(st_v_d_i_r)+length(st_v_d_i_r_incor)), length(wa_v_d_i_r)/(length(wa_v_d_i_r)+length(wa_v_d_i_r_incor)), length(pe_v_d_i_r)/(length(pe_v_d_i_r)+length(pe_v_d_i_r_incor));...
             length(st_v_d_c_s)/(length(st_v_d_c_s)+length(st_v_d_c_s_incor)), length(wa_v_d_c_s)/(length(wa_v_d_c_s)+length(wa_v_d_c_s_incor)), length(pe_v_d_c_s)/(length(pe_v_d_c_s)+length(pe_v_d_c_s_incor));...
             length(st_v_d_i_s)/(length(st_v_d_i_s)+length(st_v_d_i_s_incor)), length(wa_v_d_i_s)/(length(wa_v_d_i_s)+length(wa_v_d_i_s_incor)), length(pe_v_d_i_s)/(length(pe_v_d_i_s)+length(pe_v_d_i_s_incor))];

%% overwrite eeg stim events to match the eventlist_r
for xyz = 1:length(stimidx)
    EEGnew.event(stimidx(xyz)).switch = eventlist_r(xyz).switch;
end

%% create indexstruct= indexstruct
indexstruct.st_a_s_c = st_a_s_c;
indexstruct.wa_a_s_c = wa_a_s_c;
indexstruct.pe_a_s_c = pe_a_s_c;
indexstruct.st_a_s_i = st_a_s_i;
indexstruct.wa_a_s_i = wa_a_s_i;
indexstruct.pe_a_s_i = pe_a_s_i;
indexstruct.st_a_d_c = st_a_d_c;
indexstruct.wa_a_d_c = wa_a_d_c;
indexstruct.pe_a_d_c = pe_a_d_c;
indexstruct.st_a_d_i = st_a_d_i;
indexstruct.wa_a_d_i = wa_a_d_i;
indexstruct.pe_a_d_i = pe_a_d_i;
indexstruct.st_v_s_c = st_v_s_c;
indexstruct.wa_v_s_c = wa_v_s_c;
indexstruct.pe_v_s_c = pe_v_s_c;
indexstruct.st_v_s_i = st_v_s_i;
indexstruct.wa_v_s_i = wa_v_s_i;
indexstruct.pe_v_s_i = pe_v_s_i;
indexstruct.st_v_d_c = st_v_d_c;
indexstruct.pe_v_d_c = pe_v_d_c;
indexstruct.wa_v_d_c = wa_v_d_c;
indexstruct.st_v_d_i = st_v_d_i;
indexstruct.wa_v_d_i = wa_v_d_i;
indexstruct.pe_v_d_i = pe_v_d_i;
indexstruct.st_a_d_c_r = st_a_d_c_r;
indexstruct.wa_a_d_c_r = wa_a_d_c_r;
indexstruct.pe_a_d_c_r = pe_a_d_c_r;
indexstruct.st_a_d_c_s = st_a_d_c_s;
indexstruct.wa_a_d_c_s = wa_a_d_c_s;
indexstruct.pe_a_d_c_s = pe_a_d_c_s;
indexstruct.st_a_d_i_r = st_a_d_i_r;
indexstruct.wa_a_d_i_r = wa_a_d_i_r;
indexstruct.pe_a_d_i_r = pe_a_d_i_r;
indexstruct.st_a_d_i_s = st_a_d_i_s;
indexstruct.wa_a_d_i_s = wa_a_d_i_s;
indexstruct.pe_a_d_i_s = pe_a_d_i_s;
indexstruct.wa_a_s_c_incor = wa_a_s_c_incor;
indexstruct.pe_a_s_c_incor = pe_a_s_c_incor;
indexstruct.st_a_s_c_incor = st_a_s_c_incor;
indexstruct.st_a_s_i_incor = st_a_s_i_incor;
indexstruct.wa_a_s_i_incor = wa_a_s_i_incor;
indexstruct.pe_a_s_i_incor = pe_a_s_i_incor;
indexstruct.st_a_d_c_incor = st_a_d_c_incor;
indexstruct.wa_a_d_c_incor = wa_a_d_c_incor;
indexstruct.pe_a_d_c_incor = pe_a_d_c_incor;
indexstruct.st_a_d_i_incor = st_a_d_i_incor;
indexstruct.wa_a_d_i_incor = wa_a_d_i_incor;
indexstruct.pe_a_d_i_incor = pe_a_d_i_incor;
indexstruct.st_a_d_c_r_incor = st_a_d_c_r_incor;
indexstruct.wa_a_d_c_r_incor = wa_a_d_c_r_incor;
indexstruct.pe_a_d_c_r_incor = pe_a_d_c_r_incor;
indexstruct.st_a_d_c_s_incor = st_a_d_c_s_incor;
indexstruct.wa_a_d_c_s_incor = wa_a_d_c_s_incor;
indexstruct.pe_a_d_c_s_incor = pe_a_d_c_s_incor;
indexstruct.st_a_d_i_r_incor = st_a_d_i_r_incor;
indexstruct.wa_a_d_i_r_incor = wa_a_d_i_r_incor;
indexstruct.pe_a_d_i_r_incor = pe_a_d_i_r_incor;
indexstruct.st_a_d_i_s_incor = st_a_d_i_s_incor;
indexstruct.wa_a_d_i_s_incor = wa_a_d_i_s_incor;
indexstruct.pe_a_d_i_s_incor = pe_a_d_i_s_incor;
indexstruct.st_v_d_c_r = st_v_d_c_r;
indexstruct.wa_v_d_c_r = wa_v_d_c_r;
indexstruct.pe_v_d_c_r = pe_v_d_c_r;
indexstruct.st_v_d_c_s = st_v_d_c_s;
indexstruct.wa_v_d_c_s = wa_v_d_c_s;
indexstruct.pe_v_d_c_s = pe_v_d_c_s;
indexstruct.st_v_d_i_r = st_v_d_i_r;
indexstruct.wa_v_d_i_r = wa_v_d_i_r;
indexstruct.pe_v_d_i_r = pe_v_d_i_r;
indexstruct.st_v_d_i_s = st_v_d_i_s;
indexstruct.wa_v_d_i_s = wa_v_d_i_s;
indexstruct.pe_v_d_i_s = pe_v_d_i_s;
indexstruct.st_v_s_c_incor = st_v_s_c_incor;
indexstruct.wa_v_s_c_incor = wa_v_s_c_incor;
indexstruct.pe_v_s_c_incor = pe_v_s_c_incor;
indexstruct.st_v_s_i_incor = st_v_s_i_incor;
indexstruct.wa_v_s_i_incor = wa_v_s_i_incor;
indexstruct.pe_v_s_i_incor = pe_v_s_i_incor;
indexstruct.st_v_d_c_incor = st_v_d_c_incor;
indexstruct.wa_v_d_c_incor = wa_v_d_c_incor;
indexstruct.pe_v_d_c_incor = pe_v_d_c_incor;
indexstruct.st_v_d_i_incor = st_v_d_i_incor;
indexstruct.wa_v_d_i_incor = wa_v_d_i_incor;
indexstruct.pe_v_d_i_incor = pe_v_d_i_incor;
indexstruct.st_v_d_c_r_incor = st_v_d_c_r_incor;
indexstruct.wa_v_d_c_r_incor = wa_v_d_c_r_incor;
indexstruct.pe_v_d_c_r_incor = pe_v_d_c_r_incor;
indexstruct.st_v_d_c_s_incor = st_v_d_c_s_incor;
indexstruct.wa_v_d_c_s_incor = wa_v_d_c_s_incor;
indexstruct.pe_v_d_c_s_incor = pe_v_d_c_s_incor;
indexstruct.st_v_d_i_r_incor = st_v_d_i_r_incor;
indexstruct.wa_v_d_i_r_incor = wa_v_d_i_r_incor;
indexstruct.pe_v_d_i_r_incor = pe_v_d_i_r_incor;
indexstruct.st_v_d_i_s_incor = st_v_d_i_s_incor;
indexstruct.wa_v_d_i_s_incor = wa_v_d_i_s_incor;
indexstruct.pe_v_d_i_s_incor = pe_v_d_i_s_incor;

trimat_a = [length(st_a_s_c),length(st_a_s_c)+length(st_a_s_c_incor), length(wa_a_s_c), length(wa_a_s_c)+length(wa_a_s_c_incor), length(pe_a_s_c), length(pe_a_s_c)+length(pe_a_s_c_incor);...
             length(st_a_s_i), length(st_a_s_i)+length(st_a_s_i_incor), length(wa_a_s_i), length(wa_a_s_i)+length(wa_a_s_i_incor), length(pe_a_s_i), length(pe_a_s_i)+length(pe_a_s_i_incor);...
             length(st_a_d_c), length(st_a_d_c)+length(st_a_d_c_incor), length(wa_a_d_c), length(wa_a_d_c)+length(wa_a_d_c_incor), length(pe_a_d_c), length(pe_a_d_c)+length(pe_a_d_c_incor);...
             length(st_a_d_i), length(st_a_d_i)+length(st_a_d_i_incor), length(wa_a_d_i), length(wa_a_d_i)+length(wa_a_d_i_incor), length(pe_a_d_i), length(pe_a_d_i)+length(pe_a_d_i_incor);...
             length(st_a_d_c_r), length(st_a_d_c_r)+length(st_a_d_c_r_incor), length(wa_a_d_c_r), length(wa_a_d_c_r)+length(wa_a_d_c_r_incor), length(pe_a_d_c_r), length(pe_a_d_c_r)+length(pe_a_d_c_r_incor);...
             length(st_a_d_i_r), length(st_a_d_i_r)+length(st_a_d_i_r_incor), length(wa_a_d_i_r), length(wa_a_d_i_r)+length(wa_a_d_i_r_incor), length(pe_a_d_i_r), length(pe_a_d_i_r)+length(pe_a_d_i_r_incor);...
             length(st_a_d_c_s), length(st_a_d_c_s)+length(st_a_d_c_s_incor), length(wa_a_d_c_s), length(wa_a_d_c_s)+length(wa_a_d_c_s_incor), length(pe_a_d_c_s), length(pe_a_d_c_s)+length(pe_a_d_c_s_incor);...
             length(st_a_d_i_s), length(st_a_d_i_s)+length(st_a_d_i_s_incor), length(wa_a_d_i_s), length(wa_a_d_i_s)+length(wa_a_d_i_s_incor), length(pe_a_d_i_s), length(pe_a_d_i_s)+length(pe_a_d_i_s_incor)];

trimat_v = [length(st_v_s_c), length(st_v_s_c)+length(st_v_s_c_incor), length(wa_v_s_c), length(wa_v_s_c)+length(wa_v_s_c_incor), length(pe_v_s_c), length(pe_v_s_c)+length(pe_v_s_c_incor);...
             length(st_v_s_i), length(st_v_s_i)+length(st_v_s_i_incor), length(wa_v_s_i), length(wa_v_s_i)+length(wa_v_s_i_incor), length(pe_v_s_i), length(pe_v_s_i)+length(pe_v_s_i_incor);...
             length(st_v_d_c), length(st_v_d_c)+length(st_v_d_c_incor), length(wa_v_d_c), length(wa_v_d_c)+length(wa_v_d_c_incor), length(pe_v_d_c), length(pe_v_d_c)+length(pe_v_d_c_incor);...
             length(st_v_d_i), length(st_v_d_i)+length(st_v_d_i_incor), length(wa_v_d_i), length(wa_v_d_i)+length(wa_v_d_i_incor), length(pe_v_d_i), length(pe_v_d_i)+length(pe_v_d_i_incor);...
             length(st_v_d_c_r), length(st_v_d_c_r)+length(st_v_d_c_r_incor), length(wa_v_d_c_r), length(wa_v_d_c_r)+length(wa_v_d_c_r_incor), length(pe_v_d_c_r), length(pe_v_d_c_r)+length(pe_v_d_c_r_incor);...
             length(st_v_d_i_r), length(st_v_d_i_r)+length(st_v_d_i_r_incor), length(wa_v_d_i_r), length(wa_v_d_i_r)+length(wa_v_d_i_r_incor), length(pe_v_d_i_r), length(pe_v_d_i_r)+length(pe_v_d_i_r_incor);...
             length(st_v_d_c_s), length(st_v_d_c_s)+length(st_v_d_c_s_incor), length(wa_v_d_c_s), length(wa_v_d_c_s)+length(wa_v_d_c_s_incor), length(pe_v_d_c_s), length(pe_v_d_c_s)+length(pe_v_d_c_s_incor);...
             length(st_v_d_i_s), length(st_v_d_i_s)+length(st_v_d_i_s_incor), length(wa_v_d_i_s), length(wa_v_d_i_s)+length(wa_v_d_i_s_incor), length(pe_v_d_i_s), length(pe_v_d_i_s)+length(pe_v_d_i_s_incor)];

% write eventlist back to EEG struct
EEGold.event = eventlist;
EEG.event = eventlist_r;
EEGnew = eeg_checkset(EEG,'eventconsistency');
end % function
