function d = findContactDistance(chmap,chA,chB)



mapii = chmap == chA;
mapjj = chmap == chB;
       
[Xii, Yii] = find(mapii);
[Xjj, Yjj]  = find(mapjj);
       
d = pdist2([Xii, Yii], [Xjj, Yjj], 'euclidean');
       







return