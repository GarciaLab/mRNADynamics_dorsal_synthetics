function box = makeBox(a, b)

    box = [];

    box(1, :) = [a(1), a(2), a(3)];
    box(2, :) = [b(1), b(2), b(3)];
    box(3, :) = [a(1), b(2), a(3)];
    box(4, :) = [a(1), a(2), b(3)];
    box(5, :) = [b(1), a(2),a(3)];
    box(6, :) = [b(1), b(2),a(3)];
    box(7, :) = [a(1), b(2),b(3)];
    box(8, :) = [b(1), a(2),b(3)];



end
