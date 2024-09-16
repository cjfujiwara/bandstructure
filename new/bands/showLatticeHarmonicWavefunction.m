function [outputArg1,outputArg2] = showLatticeHarmonicWavefunction(npt,harmonic_input,opts)

if ~isfield(opts,'Indeces')
    opts.Indeces='auto';
end


% X= input.NumSites-1)

%%

if isequal(opts.Indeces,'auto')
    opts.Indeces = [];
    for kk=1:length(harmonic_input.NumBands)
        opts.Indeces(end+1)=find(harmonic_input.BandProjection(:,kk)==1,1)
    end

    % for kk=1:length(input.NumBands)-1
    %         i1 = find([input.BandProjection(:,kk)==1],1);
    %         i2 = find(input.BandProjection(:,kk+1)==1,1);
    % 
    % 
    % 
    % end
end


%%

f = figure;
f.Color='w';
f.Position=[700 100 600 400];
co =get(gca,'colororder');

% harmonic_output_H

vv = harmonic_input.EigenVectors;
pp = harmonic_input.PositionVector;
bb = harmonic_input.BandProjection;


ww = npt.Wannier_X;
wx = npt.X_extended;

xL=[-50 50];

for kk=1:length(opts.Indeces)
    %plot(harmonic_input.PositionVector,harmonic_input.EigenVectors(:,opts.Indeces(kk)),'.-');
    v = vv(:,opts.Indeces(kk));
    b = bb(opts.Indeces(kk),:);

    nband = find(round(b)==1,1);
    psi=latticeharmonic_make_position_wavefunction(v,pp,nband,ww,wx);
    

    % plot(wx,real(psi)/abs(max(psi))+opts.Indeces(kk),'-')
        plot(wx,real(psi)/max(abs(psi))+kk,'-')
        text(xL(1),kk,num2str(opts.Indeces(kk)))
        % plot(wx,imag(psi)/abs(max(psi))+kk,'-')

    % keyboard

    hold on
end
    xlim([-50 50])

end

