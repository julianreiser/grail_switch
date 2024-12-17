function[erp, erp_time, ersp, erspdb, zersp, itpc, tf_time, tf_frqs] = vs_geterpnerspnitpc(EEG, cndidx, varargin)

	%
	%
	% WHAT THIS FUNCTION DOES:
	% ------------------------
	%
    % This function takes eeg data and returns erp, db-baseline normalized ersp and itpc for each channel, as well
    % as the corresponding time and frequency vectors used for time-frequency decomposition.
	%
	% USAGE: [erp, erp_time, ersp, erspdb, itpc, tf_time, tf_frqs] = dat3d2erpnerspnitpc(EEG, {[1,2,3,4,5], [6,7,8,9,17]});
	%
	%
	% INPUT ARGUMENTS:
	% ----------------
	%
	% needed:
	% EEG              : The EEG struct
    % cndidx           : Trial indices for each condition
    %
	% optional:
	% 'blerp'          : Latencies to use for the erp baseline. Default is [-100 0].
    % 'blersp'         : Latencies to use for the ersp baseline. Default is [-500 -200].
    % 'zblersp'        : Single trial baseline latencies. Default is [-500 -200].
	% 'n_frq'          : Number of frequencies to use in tf-decomposition. Default is 20.
	% 'frqrange'       : Frequency range for tf-decomposition. Default is [3, 30].
    % 'fwhmrrange'     : Range of 'full-width-at-half-maximum' in the time domain (i.e. ms-values)
    %                    for wavelet construction. Default is [360, 100].
    % 'BL_switch'      : 1 for baseline over all conditions, 2 for condition-specific baselines
	%
	%
	% OUTPUT ARGUMENTS:
	% -----------------
	%
	% erp              : The erp matrix as channel x erp
	% erp_time         : The time vector for erp data.
    % ersp             : The ersp matrix as channel x frequency x time
    % erspdb           : The ersp matrix as channel x frequency x time decibel normalized
    % zersp            : The ersp matrix as channel x frequency x time single trial baselined (z-scores)
	% itpc             : The itpc matrix as channel x frequency x time
	% tf_time          : The time vector for time-frequency data. Differs from erp_time as edge artifacts are cropped.
	% tf_frqs          : The frequency vector for time-frequency data.
	%
	% ------------------------
	% Stefan Arnau, 27.03.2019
	% ------------------------

	% Check if any input args
	if nargin < 2
        error('Not enough input arguments... :)');
        return;
	end

	% Init input parser
	p = inputParser;

	% Set Defaults
	default_blerp = [-200, 0];
    default_blersp = [-500, -200];
    default_zblersp = [-500, -200];
	default_frqrange = [3, 30];
	default_n_frq = 28;
    default_fwhmtrange = [360, 100];
    default_resplock = 0;
    default_resplock_timewin = [-1,1];
    default_BL_switch = 2;

	% Parse inputs and set defaults
	p.FunctionName  = mfilename;
	p.CaseSensitive = false;
	p.addRequired('EEG', @isstruct);
    p.addRequired('cndidx', @iscell);
	p.addParamValue('blerp', default_blerp, @isnumeric);
    p.addParamValue('blersp', default_blersp, @isnumeric);
    p.addParamValue('zblersp', default_zblersp, @isnumeric);
	p.addParamValue('n_frq', default_n_frq, @isnumeric);
	p.addParamValue('frqrange', default_frqrange, @isnumeric);
    p.addParamValue('fwhmtrange', default_fwhmtrange, @isnumeric);
    p.addParamValue('resplock', default_resplock, @isnumeric);
    p.addParamValue('resplock_timewin', default_resplock_timewin, @isnumeric);
    p.addParamValue('BL_switch', default_BL_switch, @isnumeric);
    parse(p, EEG, cndidx, varargin{:});

    % Get times and prune
    DERP = p.Results.EEG;
    DTF = p.Results.EEG;
    pruneidx = dsearchn(DTF.times', [DTF.times(1) + 300, DTF.times(end) - 300]');
    tf_time = DTF.times(pruneidx(1) : pruneidx(2));
    erp_time = DERP.times;
    resplock = p.Results.resplock;
    resplock_timewin = p.Results.resplock_timewin;
    resplock_timewin_new = [resplock_timewin(1)*DERP.srate,resplock_timewin(2)*DERP.srate];

    % Resmats
    if resplock == 0
        erp = zeros(length(p.Results.cndidx), size(DERP.data, 1), size(DERP.data, 2));
    elseif resplock == 1
        erp = zeros(length(p.Results.cndidx), size(DERP.data, 1), length(resplock_timewin_new(1):1000/DERP.srate:resplock_timewin_new(2)));
    end % if resplock
    ersp = zeros(length(p.Results.cndidx), size(DTF, 1), p.Results.n_frq, length(tf_time));
    zersp = zeros(length(p.Results.cndidx), size(DTF, 1), p.Results.n_frq, length(tf_time));
    itpc = zeros(length(p.Results.cndidx), size(DTF, 1), p.Results.n_frq, length(tf_time));
	olmat = zeros(length(p.Results.cndidx), size(DTF, 1), p.Results.n_frq, length(tf_time));

    % Count trials in conditions and detect minimum n
    n_trials = [];
    for cnd = 1 : length(p.Results.cndidx)
        n_trials(cnd) = length(p.Results.cndidx{cnd});
    end
    mintri = min(n_trials);

    old_tricount = 1;

    % for response-locking we need to shift the ERP to the response time + the time before stim presentation
    if resplock == 1
        stim_offset = length(find(DERP.times < 0));
        timeshiftvec = ([DERP.event.rt] / (1000/DERP.srate)) + stim_offset;
    end

    % for baselining find baseline timeidx
    blerpidx = find(DERP.times >= p.Results.blerp(1) & DERP.times <= p.Results.blerp(2));

    % Iterate conditions
    for cnd = 1 : length(p.Results.cndidx)
        % Calc erps
        d = DERP.data(:, :, p.Results.cndidx{cnd});
        
        if resplock == 1
            cnd_timeshiftvec = [];
            d_timeshift = [];
            cnd_timeshiftvec = timeshiftvec(p.Results.cndidx{cnd});
            if isempty(p.Results.cndidx{cnd})
                d_timeshift = DERP.data(:,1:501,p.Results.cndidx{cnd});
            else
                for tri = 1:size(d,3)
                    d_timeshift = DERP.data(:, cnd_timeshiftvec(tri)+(resplock_timewin_new(1)*DERP.srate/1000) : cnd_timeshiftvec(tri)+(resplock_timewin_new(2)*DERP.srate/1000), p.Results.cndidx{cnd});
                    %d_timeshift = d(:,cnd_timeshiftvec(tri)+(resplock_timewin_new(1)*DERP.srate/1000) : cnd_timeshiftvec(tri)+(resplock_timewin_new(2)*DERP.srate/1000),:);
                end % for tri
            end % if isempty
        end % resplock

        for c = 1 : size(d, 1)
            fprintf('\ncondition %i/%i | Calculating ERP chan %i/%i...', cnd, length(p.Results.cndidx), c, size(d, 1));
            blerpval = [];
            if p.Results.BL_switch == 1
                blerpval = mean(squeeze(DERP.data(c, blerpidx(1):blerpidx(end), :)),[1,2]);
            elseif p.Results.BL_switch == 2
                blerpval = mean(squeeze(d(c, blerpidx(1):blerpidx(end), :)),[1,2]);
            end
            if resplock == 0
                erp(cnd, c, :) = bsxfun(@minus,mean(squeeze(d(c, :, :))', 1),blerpval);
            elseif resplock == 1
                erp(cnd, c, :) = bsxfun(@minus,mean(squeeze(d_timeshift(c, :, :))', 1),blerpval);
                erp_time = resplock_timewin_new(1):1000/DERP.srate:resplock_timewin_new(2);
            end % if resplock

            fprintf('  done');
        end

        % Set some variables
        d = DTF.data(:, :, p.Results.cndidx{cnd});
        wtime = -2 : 1 / DTF.srate : 2;
        halfw = (length(wtime) - 1) / 2;
        nconv = size(d, 2) * size(d, 3) + length(wtime) - 1;
        tf_frqs = logspace(log10(p.Results.frqrange(1)), log10(p.Results.frqrange(2)), p.Results.n_frq);
		fwhms = logspace(log10(p.Results.fwhmtrange(1)),log10(p.Results.fwhmtrange(2)),p.Results.n_frq);
        cmwX = zeros(p.Results.n_frq, nconv);

        % Build wvlts
        fprintf('\ncondition %i/%i | Building wavelets for tf-analysis...', cnd, length(p.Results.cndidx));
        for frq = 1 : p.Results.n_frq
            cmw = exp(2 * 1i * pi * tf_frqs(frq) .* wtime) .* exp((-4 * log(2) * wtime.^2) ./ (fwhms(frq) / 1000)^2); % Build
            cmw = fft(cmw ./ max(cmw), nconv); % Normalize in time domain and fft
            cmwX(frq, :) = cmw ./ max(cmw); % Normalize in frq domain
        end
        fprintf('  done');

        % Get tf trial mesures
        for c = 1 : size(d, 1)
            fprintf('\ncondition %i/%i | Calculating time-frequency decomposition chan %i/%i...', cnd, length(p.Results.cndidx), c, size(p.Results.EEG.data, 1));

            % Convolute
            pow = zeros(p.Results.n_frq, size(d, 2), size(d, 3));
            pha = zeros(p.Results.n_frq, size(d, 2), size(d, 3));
            for frq = 1 : p.Results.n_frq
                dX = fft(reshape(double(squeeze(d(c, :, :))), 1, []), nconv);
                as = ifft(cmwX(frq, :) .* dX);
                as = as(halfw + 1 : end - halfw);
                as = reshape(as, size(d, 2), size(d, 3));
                pow(frq, :, :) = abs(as) .^ 2;
                pha(frq, :, :) = angle(as);
            end

            % Cut edge artifacts
            pow = pow(:, pruneidx(1) : pruneidx(2), :);
            pha = pha(:, pruneidx(1) : pruneidx(2), :);

            % Apply single trial baseline
            blidx = dsearchn(tf_time', [p.Results.zblersp(1), p.Results.zblersp(2)]');
            zbldat = zeros(size(pow));
            for t = 1 : size(pow, 3)
                d_trial = squeeze(pow(:, :, t)); % Get trial tfmat
				for freqs = 1:size(d_trial,1)
					d_trial(freqs,:) = bsxfun(@minus,d_trial(freqs,:),mean(d_trial(freqs,:)));
					d_trial(freqs,:) = bsxfun(@rdivide, d_trial(freqs,:),std(d_trial(freqs,:)));
				end
				zbldat(:, :, t) = d_trial;
            end

            % Average
            ersp(cnd, c, :, :) = squeeze(mean(pow, 3));
            zersp(cnd, c, :, :) = squeeze(mean(zbldat,3));
            itpc(cnd, c, :, :) = abs(squeeze(mean(exp(1i*pha), 3)));

            fprintf('  done');
        end
    end % End cnd loop

    % Apply db baseline
    erspdb = zeros(size(ersp));
	zerspdb = zeros(size(zersp));
    erspga = squeeze(mean(ersp, 1));
	zerspga = squeeze(mean(zersp, 1));
    blidx = dsearchn(tf_time', [p.Results.blersp(1), p.Results.blersp(2)]');
    
    for cnd = 1 : length(p.Results.cndidx)
        for c = 1 : size(d, 1)
            
            if p.Results.BL_switch == 1
                blvals = mean(squeeze(erspga(c, :, blidx(1) : blidx(2))), 2);
			    zblvals = mean(squeeze(zerspga(c, :, blidx(1) : blidx(2))), 2);
            elseif p.Results.BL_switch == 2
                blvals = mean(squeeze(ersp(cnd,c, :, blidx(1) : blidx(2))), 2);
			    zblvals = mean(squeeze(zersp(cnd,c, :, blidx(1) : blidx(2))), 2);
            end
            
            tmp = squeeze(ersp(cnd, c, :, :));
			ztmp = squeeze(zersp(cnd, c, :, :));
            tmp = 10 * log10(bsxfun(@rdivide, tmp, blvals));
			ztmp = 10 * log10(bsxfun(@rdivide, ztmp, zblvals));
            erspdb(cnd, c, :, :) = tmp;
			zerspdb(cnd, c, :, :) = ztmp;
        end
    end
	fprintf('\n');

end
