function showBandProjections(input)
ff=figure;
ff.Color='w';
ff.Position=[50 50 800 300];
clf
co=get(gca,'colororder');

subplot(141);
plot(input.BandProjection(:,1),'.','color',co(1,:))
xlim([1 size(input.BandProjection,1)])
xlabel('eigen index')
ylabel('s projection')

subplot(142);
plot(input.BandProjection(:,2),'.','color',co(2,:))
xlabel('eigen index')
ylabel('p projection')
xlim([1 size(input.BandProjection,1)])

subplot(143);
plot(input.BandProjection(:,3),'.','color',co(3,:))
xlabel('eigen index')
ylabel('d projection')
xlim([1 size(input.BandProjection,1)])

subplot(144);
plot(input.BandProjection(:,4),'.','color',co(4,:))
xlabel('eigen index')
ylabel('f projection')
xlim([1 size(input.BandProjection,1)])

end

