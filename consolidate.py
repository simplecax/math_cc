import os
import argparse

def generate_tree(dir_path, extensions, abs_output_file_to_exclude, dirs_to_ignore, prefix=""):
    """
    Recursively generates a string representation of the file tree,
    only including specified file types and directories.
    """
    tree_lines = []
    try:
        # Get items in the directory, sorted for consistent order
        items = sorted(os.listdir(dir_path))
    except FileNotFoundError:
        return [f"ERROR: Directory not found at {dir_path}"]
    except PermissionError:
        return [f"ERROR: Permission denied to access {dir_path}"]

    # Filter out items that are not directories or don't have the right extension
    filtered_items = []
    for item in items:
        path = os.path.join(dir_path, item)

        # MODIFIED: Skip any directory present in the ignore list.
        if item in dirs_to_ignore and os.path.isdir(path):
            continue

        # Exclude the output file itself from the tree
        if os.path.abspath(path) == abs_output_file_to_exclude:
            continue

        if os.path.isdir(path):
            # Check if the directory or its subdirectories contain any valid files
            # (that are not the output file) to decide if it should be in the tree.
            should_include = False
            for sub_dirpath, sub_dirnames, sub_filenames in os.walk(path):
                # MODIFIED: Also prune ignored directories from this check for efficiency
                sub_dirnames[:] = [d for d in sub_dirnames if d not in dirs_to_ignore]
                for sub_filename in sub_filenames:
                    if any(sub_filename.endswith(ext) for ext in extensions):
                        sub_filepath = os.path.abspath(os.path.join(sub_dirpath, sub_filename))
                        if sub_filepath != abs_output_file_to_exclude:
                            should_include = True
                            break
                if should_include:
                    break
            
            if should_include:
                 filtered_items.append(item)
        elif any(item.endswith(ext) for ext in extensions):
            filtered_items.append(item)

    # Build the tree structure lines
    for i, item in enumerate(filtered_items):
        path = os.path.join(dir_path, item)
        is_last = i == (len(filtered_items) - 1)
        
        # Add the current item to the tree
        tree_lines.append(prefix + ("└── " if is_last else "├── ") + item)
        
        # If it's a directory, recurse into it
        if os.path.isdir(path):
            new_prefix = prefix + ("    " if is_last else "│   ")
            # MODIFIED: Pass the ignore list in the recursive call
            tree_lines.extend(generate_tree(path, extensions, abs_output_file_to_exclude, dirs_to_ignore, new_prefix))
            
    return tree_lines

def consolidate_project(root_dir, output_file, extensions):
    """
    Generates a file tree and consolidates the content of specified files
    into a single output text file.
    """
    # MODIFIED: Define the list of directories to ignore here for easy access.
    dirs_to_ignore = ['extern', 'build']

    if not os.path.isdir(root_dir):
        print(f"Error: The specified root directory '{root_dir}' does not exist.")
        return

    print(f"Starting consolidation for project at: {root_dir}")
    print(f"Ignoring directories: {', '.join(dirs_to_ignore)}")
    print(f"Output will be saved to: {output_file}")
    
    # Get the absolute path of the output file to avoid including it
    abs_output_file = os.path.abspath(output_file)
    
    try:
        with open(output_file, 'w', encoding='utf-8', errors='ignore') as outfile:
            # --- 1. Write the project root directory ---
            outfile.write(f"Project Root: {os.path.abspath(root_dir)}\n\n")

            # --- 2. Generate and write the file tree ---
            outfile.write("="*80 + "\n")
            outfile.write(" " * 35 + "FILE TREE\n")
            outfile.write("="*80 + "\n\n")
            # MODIFIED: Pass the ignore list to the tree generator
            tree_lines = generate_tree(root_dir, extensions, abs_output_file, dirs_to_ignore)
            for line in tree_lines:
                outfile.write(line + "\n")
            outfile.write("\n\n")

            # --- 3. Write all file contents ---
            outfile.write("="*80 + "\n")
            outfile.write(" " * 33 + "FILE CONTENTS\n")
            outfile.write("="*80 + "\n\n")

            # Walk through the directory tree
            for dirpath, dirnames, filenames in os.walk(root_dir):
                # MODIFIED: Prune ignored directories from traversal using the list.
                # This modifies the list in-place, preventing os.walk from visiting them.
                dirnames[:] = [d for d in dirnames if d not in dirs_to_ignore]

                # Sort filenames for a consistent order
                for filename in sorted(filenames):
                    if any(filename.endswith(ext) for ext in extensions):
                        file_path = os.path.join(dirpath, filename)

                        # Exclude the output file itself
                        if os.path.abspath(file_path) == abs_output_file:
                            continue
                        
                        relative_path = os.path.relpath(file_path, root_dir)
                        
                        # Write a header for each file
                        outfile.write("-" * 25 + f" FILE: {relative_path} " + "-" * 25 + "\n\n")
                        
                        try:
                            # Read file content and write it to the output
                            with open(file_path, 'r', encoding='utf-8', errors='ignore') as infile:
                                outfile.write(infile.read())
                            outfile.write("\n\n")
                        except Exception as e:
                            outfile.write(f"!!! Could not read file: {e} !!!\n\n")
        
        print("Consolidation complete!")

    except IOError as e:
        print(f"Error: Could not write to output file '{output_file}'. Reason: {e}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Consolidate project files into a single text file.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    
    parser.add_argument(
        "root_dir", 
        help="The root directory of the project you want to consolidate."
    )
    
    parser.add_argument(
        "output_file", 
        nargs='?', 
        default="consolidated_project.txt",
        help="The name of the output file. (default: consolidated_project.txt)"
    )
    
    parser.add_argument(
        "--ext",
        nargs='+',
        default=['.h', '.hpp', '.c', '.cpp', '.cxx', '.txt', '.md'],
        help="List of file extensions to include. (default: .h .hpp .c .cpp .cxx .txt .md)"
    )
    
    args = parser.parse_args()
    
    consolidate_project(args.root_dir, args.output_file, args.ext)