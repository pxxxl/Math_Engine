#pragma once
#pragma once

#include"Tools.hpp"
#include"Matrix.hpp"

namespace me {

	using variable = std::string;
	using coefficient = double;
	using term = struct { coefficient coe; variable var; };
	using polynomial = std::vector<term>;

	using equation = struct { polynomial left; polynomial right; };
	using solution = std::unordered_map<variable, polynomial>;

	using solve_status = int;

	constexpr solve_status UNSOLVED = 0;
	constexpr solve_status ERROR = 1;
	constexpr solve_status SOLVED = 2;

	class linearEquations
	{
	public:
		linearEquations(const linearEquations&) = delete;
		linearEquations(const linearEquations&&) = delete;
		linearEquations() = default;

		void inputEquation(equation a);

		void setFreeVariable(variable free_varieble);

		void solve();

		solution getSolution();

		solve_status getSolveStatus() { return status; }

		void clear();

	private:

		unsigned equation_count = 0;
		unsigned variable_count = 0;
		unsigned free_variable_count = 0;
		std::vector<variable> free_variables;
		std::vector<variable> constraint_variables;

		std::unordered_map<variable, unsigned> serial_map;
		std::vector<variable> anti_serial_map;

		std::vector<equation> equations;

		matrix<double>* solution_matrix;
		solution solution_map;

		solve_status status = UNSOLVED;
	};


	inline void linearEquations::inputEquation(equation a)
	{
		equations.push_back(a);
		equation_count++;
		auto& left_poly = a.left;
		auto& right_poly = a.right;
		auto notFindVariable = [this](std::string sample)->bool {
			for (auto& var : free_variables) {
				if (var == sample) {
					return false;
				}
			}
			for (auto& var : constraint_variables) {
				if (var == sample) {
					return false;
				}
			}
			return true;
		};
		auto insertVariebles = [this, notFindVariable](polynomial& poly) {
			for (auto& item : poly) {
				if (notFindVariable(item.var)) {
					constraint_variables.push_back(item.var);
					variable_count++;
				}
			}
		};
		insertVariebles(left_poly);
		insertVariebles(right_poly);
		return;
	}

	inline void linearEquations::setFreeVariable(variable free_varieble)
	{
		//find : iterator return
		auto findVariableFromConstrait = [this](variable sample)->decltype(constraint_variables.begin()) {
			auto be = constraint_variables.begin();
			auto en = constraint_variables.end();
			for (; be != en; be++) {
				if (*be == sample) {
					return be;
				}
			}
			return en;
		};
		//find : false return
		auto findVariableFromFree = [this](std::string sample)->int {
			for (auto& item : free_variables) {
				if (item == sample) {
					return false;
				}
			}
			return true;
		};
		decltype(free_variables.begin()) target;
		for (auto item : free_variables) {
			target = findVariableFromConstrait(item);
			if (target == constraint_variables.end()) {
				//not find 
				throw std::runtime_error("From eci::linearEquations:setFreeVariable : free varieble " + item + " set failed, no such varieble has stored in constraint variables");
			}
			else {
				if (!findVariableFromFree(item)) {
					//find in free variable
					throw std::runtime_error("From eci::linearEquations:setFreeVariable : free varieble " + item + " set failed, varieble has already stored in free variables");
				}
				constraint_variables.erase(target);
				free_variables.push_back(item);
				variable_count;
				free_variable_count++;
			}
		}
		return;
	}

	inline void linearEquations::solve()
	{
		unsigned num = 0;
		for (auto item : constraint_variables) {
			anti_serial_map.push_back(item);
			serial_map.insert(std::pair(item, num));
			num++;
		}
		for (int count = free_variables.size() - 1; count >= 0; count--) {
			anti_serial_map.push_back(free_variables[count]);
			serial_map.insert(std::pair(free_variables[count], num));
			num++;
		}
		anti_serial_map.push_back(std::string(""));
		serial_map.insert(std::pair(std::string(""), num));

		//construct matrix;
		solution_matrix = new matrix<double>(equation_count, variable_count + 1);
		unsigned cur_column = 0;
		unsigned line = 0;
		for (auto& equa : equations) {
			for (auto& terms : equa.left) {
				cur_column = serial_map[terms.var];
				(*solution_matrix).reset(line, cur_column, (*solution_matrix).get(line, cur_column) + terms.coe);
			}
			for (auto& terms : equa.right) {
				cur_column = serial_map[terms.var];
				(*solution_matrix).reset(line, cur_column, (*solution_matrix).get(line, cur_column) + terms.coe);
			}
			line++;
		}

		//solve the matrix
		*solution_matrix = (*solution_matrix).std_line_form();

		//check the result
		unsigned mat_line_max = equation_count;
		unsigned mat_column_max = variable_count + 1;
		auto is_solutionless = [&](matrix<double>* mat)->bool {
			bool flag = true;
			unsigned right_column = mat_column_max - 1;
			for (unsigned line = 0; line < mat_line_max; line++) {
				auto co = mat->get(line, right_column);
				if (is_zero(co)) {
					continue;
				}
				for (unsigned left_column = 0; left_column < right_column; left_column++) {
					co = mat->get(line, left_column);
					if (is_zero(co)) {
						continue;
					}
					flag = false;
				}
			}
			return flag;
		};
		if (is_solutionless(solution_matrix)) {
			status = ERROR;
			return;
		}
		else {
			solution_map.clear();
			unsigned cur_line = 0;
			unsigned cur_column = 0;
			while (cur_line < mat_line_max && cur_column < mat_column_max - 1) {
				if ((*solution_matrix).get(cur_line, cur_column) == 0) {
					cur_column++;
				}
				else {
					auto x_name = anti_serial_map[cur_column];
					polynomial x_solve;
					for (unsigned i = cur_column + 1; i < mat_column_max; i++) {
						auto co = solution_matrix->get(i, cur_line);
						if (!is_zero(co)) {
							x_solve.push_back(term{ co, anti_serial_map[i] });
						}
					}
					solution_map.insert(std::pair(x_name, x_solve));
					cur_line++;
				}
			}
			status = SOLVED;
		}
		return;
	}

	inline solution linearEquations::getSolution()
	{
		return solution_map;
	}

	inline void linearEquations::clear()
	{
		unsigned equation_count = 0;
		unsigned variable_count = 0;
		unsigned free_variable_count = 0;
		free_variables.clear();
		constraint_variables.clear();
		serial_map.clear();
		equations.clear();
		solution_matrix->clear();
		solution_map.clear();
		status = UNSOLVED;
	}

}


