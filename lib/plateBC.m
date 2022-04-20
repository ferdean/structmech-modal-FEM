function [constrNodes]=plateBC(nodes,side,floor)



nN          = size(nodes,1);
constrNodes = [ ];

if side == 'axis'
    for n = 1:nN
        if nodes(n,1) == floor
            constrNodes = [constrNodes n]; 
        end
    end
elseif side == 'base'
    for n = 1:nN
        if nodes(n,2) == floor
            constrNodes = [constrNodes n]; 
        end
    end
end

constrNodes=constrNodes.';

end

