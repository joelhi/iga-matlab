function index = findIndexJoel(vec,dof)
lenVec = length(vec);

for i = 1:lenVec
    dofVec = vec(i);
    if abs(dof - dofVec) < 0.0001
    index = i;
    end
end
end
