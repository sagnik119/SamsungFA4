load Data_SIMO.mat

bu_a_lin = [];

for i = 1:numel(Data_SIMO)
    data = Data_SIMO(i);
    bu_a_lin = [bu_a_lin; data.bu_a_lin];
end