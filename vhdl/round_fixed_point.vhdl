
library ieee;

use ieee.std_logic_1164.all;
use std.textio.all;

entity round_fixed_point_testbench is
end entity round_fixed_point_testbench;

architecture arch of round_fixed_point_testbench is

    signal CLK : std_logic;
    signal sig : std_logic_vector(7 downto 0);

    function as_string(slv : in std_logic_vector) return string is
    begin
        return string'("hallo");
    end function as_string;

begin

    process (CLK) is
        variable l : line;
    begin
        if rising_edge(CLK) then
            write(l, string'("sig: "));
            write(l, as_string(sig));
            writeline(output, l);
        end if;
    end process;

    process is
        variable l : line;
    begin
        for i in 1 to 100 loop
            write(l, string'("i: "));
            write(l, i);
            writeline(output, l);
            CLK <= '1';
            wait for 10 ns;
            CLK <= '0';
            wait for 10 ns;
        end loop;
        wait;
    end process;

end architecture arch;
