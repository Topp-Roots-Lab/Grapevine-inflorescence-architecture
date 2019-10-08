function adj=edgeL2adj(el)
nodes=sort(unique([el(:,1) el(:,2)])); % get all nodes, sorted
adj=zeros(numel(nodes));   % initialize adjacency matrix
% across all edges
for i=1:size(el,1);adj(find(nodes==el(i,1)),find(nodes==el(i,2)))=el(i,3);
end