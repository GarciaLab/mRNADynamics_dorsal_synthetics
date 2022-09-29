spotIndex = 3;
s = spot_struct_protein(spotIndex)
s.setID
s.nucleusID
set.ncID
s.frames
s.edge_qc_flag_vec
badEdgeFrames = s.frames(s.edge_qc_flag_vec==-1);
badSerialFrames = s.frames(s.serial_qc_flag_vec==-1);
snips = snip_data([snip_data.setID]==s.setID)

for k = 1:length(badEdgeFrames)
    row = [snips.frame]==badEdgeFrames(k);
    f = figure; tiledlayout('flow');
    nexttile;
    imagesc(imresize(snips(row).spot_protein_snips, 10));
    nexttile
    imagesc(imresize(snips(row).spot_mcp_snips, 10));
    nexttile
    imagesc(imresize(snips(row).edge_control_protein_snips, 10));
    nexttile
    imagesc(imresize(snips(row).edge_control_mcp_snips, 10));
    waitforbuttonpress;
    close(f)
end


for k = 1:length(badSerialFrames)
    row = [snips.frame]==badSerialFrames(k);
    f = figure; tiledlayout('flow');
    nexttile;
    imagesc(imresize(snips(row).spot_protein_snips, 10));
    nexttile
    imagesc(imresize(snips(row).spot_mcp_snips, 10));
    nexttile
    imagesc(imresize(snips(row).edge_control_protein_snips, 10));
    nexttile
    imagesc(imresize(snips(row).edge_control_mcp_snips, 10));
    waitforbuttonpress;
    close(f)
end