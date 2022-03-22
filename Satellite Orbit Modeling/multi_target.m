% obser = load('obs_rv.mat');
% obs_rv = obser.obs_rv; 
% obs_k = obs_rv(1,:);


multiobj = load('obs_simulation_rv.mat');
multiobj = multiobj.obs_target1;
more_eph = zeros(size(multiobj,1),241,6);
for k = 1:size(multiobj,1)
    Eph_eci = orbit_prop(multiobj(k,:));
    more_eph(k,:,:) = Eph_eci;
end

save('more_eph.mat')