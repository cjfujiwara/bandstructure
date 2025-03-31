function wfs = calculateLatticeHarmonicWavefunction(npt,harmonic_input,opts)

if ~isfield(opts,'Indeces')
    opts.Indeces='auto';
end



%%

if isequal(opts.Indeces,'auto')
    opts.Indeces = [];
    for kk=1:length(harmonic_input.NumBands)
        opts.Indeces(end+1)=find(harmonic_input.BandProjection(:,kk)==1,1)
    end
end


%%


vv = harmonic_input.EigenVectors;
pp = harmonic_input.PositionVector;
bb = harmonic_input.BandProjection;
ww = npt.Wannier_X;
wx = npt.X_extended;

wfs = zeros(numel(wx),length(opts.Indeces));

for kk=1:length(opts.Indeces)
    v = vv(:,opts.Indeces(kk));
    b = bb(opts.Indeces(kk),:);
    nband = find(round(b)==1,1);
    psi=latticeharmonic_make_position_wavefunction(v,pp,nband,ww,wx);    
    wfs(:,kk)=psi;
end




end

